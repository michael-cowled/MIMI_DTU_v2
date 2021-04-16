##########################################################
##########################################################
# Microbial Interaction Metabolite Integrator - Script 1 #
##########################################################
##########################################################
 
# COPYWRIGHT: © Macquarie University - Michael Cowled and contributors [2021]

# AUTHOR COMMENT: MIMI is a program developed as part of my PhD and was supported
# by MQ RTP Scholarship.

# PURPOSE:For use in comparing the chemical extracts of microbial cocultures to 
# their axenic monocultures.

# INPUT: Interaction_Matrix - A df containing the columns CON1, CON2	& Coculture
# is required to be made and imported by the user to compare samples.

# OUTPUT: An output file (.csv) is created matching the peaks from the coculture
# to peaks in either of the two controls, quantifying the change in peak area.
# Further manipulation is needed in Script_2 to remove doubley assigned peaks,
# identify non UVs and add in peaks belonging to the control(s).

# R packages required to be loaded in:

library(dplyr)
library(tidyr)
library(readxl)

#------------------------------------------------------------------------------#

# Functions to be pre-loaded prior to use of the main function, MIMI()

# 1.ReadExcel - Reads in the NovaC files corresponding to a particular row number 
# in the Interaction_Matrix. See MIMI_Code_Book for explanation of NovaC.

# 2.CheckUVCount - Compares the top 5 UV maxima of a matched peak in the control 
# of interest and the coculture.

# 3.ReadUV - Reads in the UV spectral data (abs vs. wavelength) for the 
# corresponding 3 samples being compared.

# 4.SubtractUV - Subtracts the UV spectrum of the control of interest from the 
# coculture, and takes the mean(Abs).

# 5.RowBinder - Binds the matched peak found in PeakMatcher to a df named cc.df

# 6.CalcRatio: Calculates the %Enhancement/Suppresion compared to control levels

# 7.PeakMatcher - PeakMatcher: Matches and verifies peaks from the control of 
# interest to peaks in the coculture, utilising a combination of retention time, 
# number of matching UV maxima (UVcheck) and the means of the subtracted UV 
# spectra (UVsubtract).

# 8.ConConsolidator - Consolidates the outcome table to match peaks from coculture 
# to a single peak in a control. Compares subtracted UV spectra to make decisions 
# based on double matching.

# 9.EffectCategoriser - Characterises the Peak Area ratio as an effect to the 
# metabolite in the coculture (induction, suppression, etc.).

#------------------------------------------------------------------------------#

# Special Terms:
# 1. CC - coculture: the duoculture of two organisms grown together.
# 2. con - control: the axenic monoculture of one of the organisms.

# Reading in the functions:

#############################################
# 1.ReadExcel
#############################################

#Reads in NovaC excel file and reformats UV column to be more usable.

#Arg: 
# excel.name is read in from Interaction_Matrix object as CC or CON name.

#Returns:
# raw.df.with.sorted.uv - a df corresponding to the excel file read in.

ReadExcel <- function(excel.name) {

    # Read in a df based on NovaC excel format
    raw.df <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", 
                                excel.name, ".xlsm"), skip = 3)
    
    # A new df is set up to capture the top 5 UV maxima for each peak
    uv.df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
                      c("UV1", "UV2", "UV3", "UV4", "UV5"))
    
    # The following extracts and sorts the top 5 UV maxima from each peak
    n <- nrow(raw.df)
    cc.peak <- 1
    while (cc.peak <= n) {
        if (!is.na(raw.df$'UV Peaks')) {
            raw.uv <-raw.df$`UV Peaks`[[cc.peak]]
            raw.uv.processed <- gsub("\\(","", raw.uv)
            raw.uv.processed <- gsub("\\)","", raw.uv.processed)
            raw.uv.processed <- gsub("< 190","", raw.uv.processed)
            raw.uv.processed <- gsub("s","", raw.uv.processed)
            raw.uv.split <- strsplit(raw.uv.processed, "\r\n")
            raw.uv.split <- as.data.frame(raw.uv.split, col.names = "UV1")
            raw.uv.split <- separate(raw.uv.split, 'UV1', c("UV1", "UV1_percent"), 
                                     sep = " ")
            raw.uv.split <- transform(raw.uv.split, 
                                      UV1_percent = as.numeric(UV1_percent))
            raw.uv.ordered <- raw.uv.split[order(raw.uv.split$UV1_percent, 
                                                 decreasing = TRUE), ]
            
            # A second uv.df is set up to capture the UV data for a single peak 
            uv.df.peak <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), 
                                   c("UV1", "UV2", "UV3", "UV4", "UV5"))
            
            # Transposes UVs into row format; more efficient method possible
            uv.df.peak[1, 1] <- raw.uv.ordered[1, 1]
            uv.df.peak[1, 2] <- raw.uv.ordered[2, 1]
            uv.df.peak[1, 3] <- raw.uv.ordered[3, 1]
            uv.df.peak[1, 4] <- raw.uv.ordered[4, 1]
            uv.df.peak[1, 5] <- raw.uv.ordered[5, 1]
            
            # Lists the UV for the associated peak into the originally setup uv.df
            uv.df <- rbind(uv.df, uv.df.peak)
            
        }   else {
            
            # Lists the UVs as nothing instead
            uv.df.peak <- c(NA, NA, NA, NA, NA)
            uv.df <- rbind(uv.df, uv.df.peak)
            
        }   
        cc.peak <- cc.peak + 1
    }
    
    # uv.df is added onto the loaded raw.df object, correlating peaks with UV
    names(uv.df) <- c("UV1", "UV2", "UV3", "UV4", "UV5")   
    uv.df <- transform(uv.df, UV1 = as.numeric(UV1))
    uv.df <- transform(uv.df, UV2 = as.numeric(UV2))
    uv.df <- transform(uv.df, UV3 = as.numeric(UV3))
    uv.df <- transform(uv.df, UV4 = as.numeric(UV4))
    uv.df <- transform(uv.df, UV5 = as.numeric(UV5))
    raw.df.with.sorted.uv <- cbind(raw.df, uv.df)
    raw.df.with.sorted.uv <- select(raw.df.with.sorted.uv, 1:4, 20:24)
    return(raw.df.with.sorted.uv)
    
}

#############################################
# 2.CheckUVCount
#############################################

# Compares the number of matching UV maxima in the con to cc.
# UV data is located in columns 5:9, hence 'con.uv.no' and 'cc.uv.no' starts at 5.

# Args:
    # con can be either the df generated for control 1 or control 2.
    # cc is the df generated for the coculture.
    # cc.peak is the coculture peak no. to be compared.
    # con.peak is the control peak no. to be compared.

# Returns:
    # uv.count - the number of matching uv maxima for the peaks being compared
    # in the control and the coculture.

CheckUVCount <- function(con, cc, cc.peak, con.peak) {
    
    uv.count <- 0
    cc.uv.no <- 5 
    con.uv.no <- 5
    
    while (cc.uv.no < 10) {
        con.uv.no <- 5
        while (con.uv.no < 10 && cc.uv.no < 10) {
            if (con.uv.no == 9)  {
                con.uv.no <- con.uv.no + 1
                cc.uv.no <- cc.uv.no + 1
            }   else if (is.na(cc[cc.peak, cc.uv.no])) {
                cc.uv.no <- cc.uv.no + 1
            }   else if (is.na(con[con.peak, con.uv.no])) {
                con.uv.no <- con.uv.no + 1
            }   else if (cc[cc.peak, cc.uv.no] <= con[con.peak, con.uv.no] + 2 && 
                         cc[cc.peak, cc.uv.no] >= con[con.peak, con.uv.no] - 2) {
                uv.count <- uv.count + 1
                con.uv.no <- con.uv.no + 1
            }   else {
                con.uv.no <- con.uv.no + 1
            }
        }
    }
    return(uv.count)
}

#############################################
# 3.ReadUV Function
#############################################

# Interprets and reads the raw UV data for every wavelength.

# Args:
    # excel.name is the name of the cc or con read in from the Interaction_Matrix.

# Returns:
    # uv.name - a df with col1 = Wavelength, col2 = Abs for peak 1, etc.

ReadUV <- function(excel.name) {
    uv.name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/",
                                 excel.name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(uv.name)
}

#############################################
# 4.SubtractUV: 
#############################################

# Subtracts whole UV spectra and calculates the mean for peaks being compared
# in the control and coculture. The resulting subtraction will provide a value
# close to ZERO if the UV spectra are close in shape.

# Args:
    # con.uv is the uv df generated for the control being compared.
    # cc.uv is the uv df generated for the coculture being compared.
    # cc.peak is the peak no. being compared in the coculture.
    # con.peak is the peak no. being compared in the control.

# Returns:
    # uv.mean - the mean of UV absorbances of the con subtracted from the cc.


SubtractUV <- function(con.uv, cc.uv, cc.peak, con.peak) {
    
    cc.uv.no <- cc.peak + 1
    con.uv.no <- con.peak + 1
    
    # 'cc.uv.no' refers to the column of abs data being compared in coculture
    # 'con.uv.no' refers to the column of abs data being compared in control
    # Note: Col2 = Peak#1, Col3 = Peak#2, etc. hence '+1'
    
    subtracted.uv <- cc.uv[, cc.uv.no] - con.uv[, con.uv.no]
    uv.mean <- mean(subtracted.uv)
    return(uv.mean)
}

#############################################
# 5.RowBinder
#############################################

# Adds a row with the appropriate matched peak values.

# Args:
    # cc.df is the core df being worked on to generate the end output file.
    # coculture is the generated df for the coculture being compared.
    # con is the generated df for the control for which the peak is matching.
    # final.count is the number of matching maxima between the matched peaks.
    # uv.mean is the average of the subtracted UV specral absorbances.
    # ratio is the computed ratio of peak areas for the matched peaks.
    # cc.peak is the peak no. being compared in the coculture.
    # con.peak is the peak no. being compared in the control.
    # check.else is a defining parameter to indicate whether the peaks match.

# Returns:
    # cc.df - an updated cc.df with a new row indicating the matched peaks.

RowBinder <- function(cc.df, coculture, con, final.count, uv.mean, ratio, 
                      cc.peak, con.peak, check.else) {

    if (check.else == TRUE) {
        cc.df <- rbind(cc.df, 
                       c(PeakNo_CC = coculture$Peak[cc.peak], 
                         RetTime_CC = coculture$RetTime[cc.peak], 
                         PeakArea_CC = coculture$Area[cc.peak], 
                         PeakNo_con = NA, RetTime_con = NA, 
                         PeakArea_con = NA, UV_Count_con = NA, 
                         Subtracted_UV_Mean_con = NA, 
                         PeakRatio_con = NA))
    }
    else {
        cc.df <- rbind(cc.df, 
                       c(PeakNo_CC = coculture$Peak[cc.peak], 
                         RetTime_CC = coculture$RetTime[cc.peak], 
                         PeakArea_CC = coculture$Area[cc.peak], 
                         PeakNo_con = con$Peak[con.peak], 
                         RetTime_con = con$RetTime[con.peak],
                         PeakArea_con = con$Area[con.peak], 
                         UV_Count_con = final.count, 
                         Subtracted_UV_Mean_con = uv.mean, 
                         PeakRatio_con = ratio))
    }
    return(cc.df)
}

#############################################
# 6.CalcRatio
#############################################

# Calculates the %Enhancement/Suppresion compared to control levels.

# Args:
    # cc is the generated df of the coculture.
    # cc.peak is the peak number of the matched peak in the coculture.
    # con is the generated df of the matching control.
    # con.peak is the peak number of the matched peak in the control.

# Returns:
    # ratio - The ratio of peak ares of the matched peaks as a percentage.

CalcRatio <- function(cc, cc.peak, con, con.peak) {
    
    ratio <- (((cc[cc.peak, 3] -con[con.peak, 3]) / 
                   con[con.peak, 3]) * 100)
    return(ratio)
    
}

#############################################
# 7.PeakMatcher
#############################################

# Finds the nearest matching peak in control based on retention time.
# Then verifies whether the matching peak is valid based on number of matching
# UV maxima and/or the subtracted UV spectra of the two peaks is close to 0.

# Args:
    # con is the generated df of the control to be compared.
    # coculture is the generated df of the coculture to be compared.
    # con.uv is the generated df of UV absorbaces for each peak in th con.
    # cc.uv is the generated df of UV absorbances for each peak in the cc.

# Returns:
    # cc.df - a generated df to catalogue the matching peaks in the cc and con.

PeakMatcher <- function(con, coculture, con.uv, cc.uv) {
    
    # Sets up a df to with the desired column names.
    cc.df <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                      c("PeakNo_CC", "RetTime_CC", "PeakArea_CC", "PeakNo_con", 
                        "RetTime_con", "PeakArea_con", "UV_Count_con", 
                        "Subtracted_UV_Mean_con", "PeakRatio_con"))
    
    n <- nrow(coculture)
    cc.peak <- 1
    
    # 'cc.peak' corresponds to the peak no. to be compared in the coculture
    # 'con.peak' finds the peak in con1 with the closest RetTime to peak 'cc.peak'
    
    while (cc.peak < n + 1) {
        con.peak <- which(abs(con$RetTime - coculture$RetTime[cc.peak]) ==
                              min(abs(con$RetTime - coculture$RetTime[cc.peak])))
        ratio <- CalcRatio(coculture, cc.peak, con, con.peak)
        final.count <- CheckUVCount(con, coculture, cc.peak, con.peak)
        uv.mean <- SubtractUV(con.uv, cc.uv, cc.peak, con.peak)
        
        # Performs both the CheckUVCount and SubtractUV functions
        
        check.else <- FALSE  # A parameter set up for RowBinder function
        
        if (coculture$RetTime[cc.peak] < (con$RetTime[con.peak] + 0.15) && 
            coculture$RetTime[cc.peak] > (con$RetTime[con.peak] - 0.15) &&
            (final.count > 2 || abs(uv.mean) < 1.5)) {
            
            # Checks if the ret.time is within 0.15 min of each other
            # Checks the uv.count is at least 3 OR uv.mean < 1.5
            # If satisfied, the peak is declared a match and added to df

            cc.df <- RowBinder(cc.df, coculture, con, final.count, uv.mean, 
                               ratio, cc.peak, con.peak, check.else)
            
        }    else if (coculture$RetTime[cc.peak] < (con$RetTime[con.peak] + 0.15) && 
                      coculture$RetTime[cc.peak] > (con$RetTime[con.peak] -0.15) &&
                      final.count > 1 && abs(uv.mean) < 2) {  
            
            # Checks if the ret.time is within 0.15 min of each other
            # Checks the uv.count is at least 2 AND uv.mean < 2
            # Both uv.count and uv.mean are less strict requirement
            # But requiring satisfaction of both adds stringency
            # If satisfied, the peak is declared a match and added to df

            cc.df <- RowBinder(cc.df, coculture, con, final.count, uv.mean, 
                               ratio, cc.peak, con.peak, check.else)
            
        }    else if (coculture$RetTime[cc.peak] < (con$RetTime[con.peak] + 0.05) && 
                      coculture$RetTime[cc.peak] > (con$RetTime[con.peak] -0.05) &&
                      final.count > 1) {  
            
            cc.df <- RowBinder(cc.df, coculture, con, final.count, uv.mean, 
                               ratio, cc.peak, con.peak, check.else)
            
        }   else  { 
            
            # The peak is declared NOT to have a match in the control and given NA
            check.else <- TRUE
            cc.df <- RowBinder(cc.df, coculture, con, final.count, uv.mean, 
                               ratio, cc.peak, con.peak, check.else)
        }
        check.else <- FALSE
        cc.peak <- cc.peak + 1
    }
    return(cc.df)
}

#############################################
# 8.ConConsolidator
#############################################

# The following code removes double assignments to a coculture peak
# from two different controls.

# An example is that perhaps peak 1 from the coculture matches
# to peak 1 in con1 and peak 2 in con2.

# This function will determine which control's matched peak matches
# closest and remove the assignment for the weakest match.

# Args:
    # cc.df is df being worked on summarising all matched peaks in cc and cons.
    # con1.df is the generated df for control 1.
    # con2.df is the generated df for control 2.

# Returns:
    # matched.peak.df - a tidied version of cc.df matching to single con peaks.

ConConsolidator <- function(cc.df, con1.df, con2.df) {
    
    row.no <- 1
    total.rows <- nrow(cc.df)
    
    # sets up a new tidier version of the df
    
    matched.peak.df <- setNames(data.frame(matrix(ncol = 7, nrow = nrow(cc.df))), 
                                c("Matched_con", "PeakNo_con", "RetTime_con", 
                                  "PeakArea_con", "UV_Count", 
                                  "Subtracted_uv.mean", "PeakRatio"))
    
    while (row.no <= total.rows + 1) {
        if (is.na(cc.df[row.no, 4]) && is.na(cc.df[row.no, 10])) {
        }   else if (!is.na(cc.df[row.no, 4]) && 
                     !is.na(cc.df[row.no, 10] && 
                     abs(cc.df[row.no, 8]) < abs(cc.df[row.no, 14]))) { 
            
            # When there are two peaks matched (!is.na for both),
            # the uv.means are compared, with the higher removed.
            
            matched.peak.df[row.no, 1] <- con1.df
            matched.peak.df[row.no, 2:7] <- cc.df[row.no, 4:9]
        }   else if (!is.na(cc.df[row.no, 4]) && 
                     !is.na(cc.df[row.no, 10] && 
                     abs(cc.df[row.no, 8]) > abs(cc.df[row.no, 14]))) { 
            matched.peak.df[row.no, 1] <- con2.df
            matched.peak.df[row.no, 2:7] <- cc.df[row.no, 10:15]
        }   else if (is.na(cc.df[row.no, 10]))   {
            
            # If no double-peak mactching but signifies a matched con1 peak
            # then con1 peak set as the matched peak
            
            matched.peak.df[row.no, 1] <- con1.df
            matched.peak.df[row.no, 2:7] <- cc.df[row.no, 4:9]
        }   else {
            
            #If no double-peak mactching or matched con1 peak
            #then con2 peak set as the matched peak
            
            matched.peak.df[row.no, 1] <- con2.df
            matched.peak.df[row.no, 2:7] <- cc.df[row.no, 10:15]
        }
        row.no <- row.no + 1
    }
    return(matched.peak.df)
}

#############################################
# 9.Metabolite Effect Characteriser
#############################################

# Characterises the ratios of peak area changes to a categorical factor.

# Args:
    # refined.cc.df is the tidied version of the df being worked on summarising 
    # all matched peaks in cc and cons. 
    # cc.name is the name of the coculture being investigated.

#Returns:
    # refined.cc.df - includes a new column identifies matched peaks as inductions,
    # suppressions, etc.

EffectCategoriser <- function(refined.cc.df, cc.name) {
    
    row.no <- 1
    total.rows <- nrow(refined.cc.df)
    
    # A df is created to list the effects corresponding to matched peaks
    
    metabolite.effect.df <- setNames(data.frame(matrix(ncol=1, nrow=total.rows)), 
                                     c("Metabolite_Effect"))
    
    while (row.no <= total.rows) {
        if (refined.cc.df[row.no, 1] == cc.name && 
            is.na(refined.cc.df[row.no, 11])) {
            metabolite.effect.df[row.no, 1] <- 6  # Induction
        }   else if (refined.cc.df[row.no, 1] == cc.name && 
                     refined.cc.df[row.no, 11] > -100
                     && refined.cc.df[row.no, 11] <= -20) {
            metabolite.effect.df[row.no, 1] <- 2  # Suppression
        }   else if (refined.cc.df[row.no, 1] == cc.name && 
                     refined.cc.df[row.no, 11] > -20
                     && refined.cc.df[row.no, 11] < 20) {
            metabolite.effect.df[row.no, 1] <- 3  # Little to No Change
        }   else if (refined.cc.df[row.no, 1] == cc.name && 
                     refined.cc.df[row.no, 11] >= 20
                     && refined.cc.df[row.no, 11] < 100) {
            metabolite.effect.df[row.no, 1] <- 4  # Enhancement
        }   else if (refined.cc.df[row.no, 1] == cc.name && 
                     refined.cc.df[row.no, 11] >= 100) {
            metabolite.effect.df[row.no, 1] <- 5  # Major Enhancement
        }
        row.no <- row.no + 1
    }
    refined.cc.df <- cbind(refined.cc.df, metabolite.effect.df)
    return(refined.cc.df)
}

#############################################
# 9.CalcPercArea
#############################################

# Calculates the %Peak Area for the peaks in coculture
# The purpose of this function is to give induced peaks a quantifiable parameter

# Args:
    # refined.cc.df is the tidied version of the df being worked on summarising 
    # all matched peaks in cc and cons. 

# Returns:
    # refined.cc.df - includes a new column with the calculated PercArea of peaks.

CalcPercArea <- function(refined.cc.df) {
    refined.cc.df <- mutate(refined.cc.df, PercArea = PeakArea_CC / 
                                sum(refined.cc.df$PeakArea_CC) * 100)
    refined.cc.df <- refined.cc.df[c(1:4, 13, 5:12)]
    return(refined.cc.df)
}


#############################################
# MIMI_Part_1
#############################################

# Microbial Interaction Metabolite Integrator.
# The main function to peak match and compare the coculture to the two controls.

# Args:
    # Interaction_Matrix is a user defined matrix contraining the column names:
        # CON1
        # CON2	
        # Coculture
    # Note: NovaC files (a refined version of the raw HPLC data) should be 
    # present for each control and coculture.

# Returns:
    # Output files derived from refined.cc.df which summarise all matched peaks
    # in each coculture and its corresponding controls.

MIMI <- function() {           
    
    matrix.total.rows <- nrow(Interaction_Matrix)
    matrix.row.no <- 1
    
    while (matrix.row.no <= matrix.total.rows) {
        
        con1.df <- as.character(Interaction_Matrix[matrix.row.no, 1])
        con2.df <- as.character(Interaction_Matrix[matrix.row.no, 2])
        cc.name <- as.character(Interaction_Matrix[matrix.row.no, 3])
        matrix.row.no <- matrix.row.no + 1
        print(cc.name)
        #Generates 3 dataframes using the ReadExcel function for the first
        #interaction to be investigated from the Interaction_Matrix
        
        con1 <- as.data.frame(ReadExcel(con1.df))
        con2 <- as.data.frame(ReadExcel(con2.df))
        coculture <- as.data.frame(ReadExcel(cc.name))
        
        # Generates 3 dataframes using the ReadUV function
        
        con1.uv <- as.data.frame(ReadUV(con1.df))
        con2.uv <- as.data.frame(ReadUV(con2.df))
        cc.uv <- as.data.frame(ReadUV(cc.name))
        
        # Runs the peak matching algorithm for con1 and con2 separately.
        # Then merges the two together into a unified df.
        
        cc.df <- PeakMatcher(con1, coculture, con1.uv, cc.uv)
        cc.df2 <- PeakMatcher(con2, coculture, con2.uv, cc.uv)
        cc.df.merged <- cbind(cc.df, cc.df2[, 4:9])
        
        # Performs a function to correct for double peak matching to a unique
        # peak to peaks from more than one con
        
        matched.peak.df <- ConConsolidator(cc.df.merged, con1.df, con2.df)
        
        # Final processing steps in creating the tidied df (refined.cc.df)
        
        matched.peak.df <- matched.peak.df[1:nrow(matched.peak.df), ]
        df <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(cc.df.merged))), 
                       "Sample_Ref")
        df[1:nrow(df), ] <- cc.name
        refined.cc.df <- cc.df[1:nrow(cc.df.merged), 1:3]
        refined.cc.df <- cbind(df, refined.cc.df)
        refined.cc.df <- cbind(refined.cc.df, matched.peak.df)
        
        # Performs a function to characterise effects onto metabolites
        # and adds this into the tidied data set
        
        refined.cc.df <- 
            EffectCategoriser(refined.cc.df, cc.name)
        names(refined.cc.df)[2:4] <- c("PeakNo_CC", "RetTime_CC", "PeakArea_CC")
        refined.cc.df <- CalcPercArea(refined.cc.df)
        
        # Rewrites the tidied dataset to a csv file
        
        write.csv(refined.cc.df, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         cc.name, ".CSV"), row.names = FALSE)
    }
}

# Load in excel file named Interaction_Matrix and use MIMI() function:
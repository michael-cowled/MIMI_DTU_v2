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

# 5.PeakMatcher - PeakMatcher: Matches and verifies peaks from the control of 
# interest to peaks in the coculture, utilising a combination of retention time, 
# number of matching UV maxima (UVcheck) and the means of the subtracted UV 
# spectra (UVsubtract).

# 6.ConConsolidator - Consolidates the outcome table to match peaks from coculture 
# to a single peak in a control. Compares subtracted UV spectra to make decisions 
# based on double matching.

# 7.EffectCategoriser - Characterises the Peak Area ratio as an effect to the 
# metabolite in the coculture (induction, suppression, etc.).

#------------------------------------------------------------------------------#

# Special Terms:
# 1. CC - coculture: the duoculture of two organisms grown together.
# 2. con - control: the axenic monoculture of one of the organisms.

# Reading in the functions:

#############################################
# 1.ReadExcel: Creates 3 data frames for the 3 samples to be compared.
#############################################

ReadExcel <- function(excel.name) {

#Reads in NovaC excel file and reformats UV column to be more usable.

#Arg: 
    #excel.name is read in from Interaction_Matrix object as CC or CON name.
    
#Returns:
    #A df corresponding to the excel file read in.
    
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
# 2.CheckUVCount: Verifies the matching UV maxima for peaks in con1/2 & coculture
#############################################

# CheckUVCount: Comparing con1 to coculture
# Inputs are con1/con2 and coculture as read from the Interaction_Matrix
# 'cc.peak' refers to the peak no. being compared in coculture
# 'con.peak' refers to the peak no. being compared in control
# UV data is located in columns 5:9, hence 'con.uv.no' and 'cc.uv.no' starts at 5.

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
# 3.ReadUV Function: Reads in the raw UV data for every wavelength
#############################################

# Input is excel.name which is read from the Interaction_Matrix object
# Generates a df with col1 = Wavelength, col2 = Abs for Peak 1, etc.

ReadUV <- function(excel.name) {
    uv.name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/",
                                 excel.name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(uv.name)
}

#############################################
# 4.SubtractUV: Subtracts whole UV spectra and calculates the mean for peaks 
# in con1/2 & coculture
#############################################

# SubtractUV: Comparing control to coculture

SubtractUV <- function(con.uv, cc.uv, cc.peak, con.peak) {
    
    # Inputs are con1.uv and cc.uv as created from the Interaction_Matrix
    # 'cc.peak' refers to the peak no. being compared in coculture
    # 'con.peak' refers to the peak no. being compared in control
    
    cc.uv.no <- cc.peak + 1
    con.uv.no <- con.peak + 1
    
    # Subtracts absorbances at each wavelength and calculates the column mean
    # 'cc.uv.no' refers to the column of abs data being compared in coculture
    # 'con.uv.no' refers to the column of abs data being compared in control
    # Note: Col2 = Peak#1, Col3 = Peak#2, etc. hence '+1'
    
    subtracted.uv <- cc.uv[, cc.uv.no] - con.uv[, con.uv.no]
    uv.mean <- mean(subtracted.uv)
    return(uv.mean)
}

#############################################
# 5.PeakMatcher: Finds and compares the nearest matching peak in control
#############################################

PeakMatcher <- function(con, coculture, con.uv, cc.uv) {
    
    # Sets up a df to with the desired column names.
    cc.df <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                      c("PeakNo_CC", "RetTime_CC", "PeakArea_CC", "PeakNo_con", 
                        "RetTime_con", "PeakArea_con", "UV_Count_con", 
                        "Subtracted_UV_Mean_con", "PeakRatio_con"))
    
    n <- nrow(coculture)
    cc.peak <- 1
    
    # 'cc.peak' corresponds to the peak no. to be compared in the coculture
    
    while (cc.peak < n + 1) {
        con.peak <- which(abs(con$RetTime - coculture$RetTime[cc.peak]) ==
                              min(abs(con$RetTime - coculture$RetTime[cc.peak])))
        
        # 'con.peak' finds the peak in con1 with the closest RetTime to peak 'cc.peak'
        
        ratio <- (((coculture[cc.peak, 3] - con[con.peak, 3]) / 
                       con[con.peak,3]) * 100) 
        
        # Computes the ratio of peak areas as a %
        
        final.count <- CheckUVCount(con, coculture, cc.peak, con.peak)
        uv.mean <- SubtractUV(con.uv, cc.uv, cc.peak, con.peak)
        
        # Performs both the CheckUVCount and SubtractUV functions
        
        if (coculture$RetTime[cc.peak] < (con$RetTime[con.peak] + 0.15) && 
            coculture$RetTime[cc.peak] > (con$RetTime[con.peak] - 0.15) &&
            (final.count > 2 || abs(uv.mean) < 1.5)) {
            
            # Checks if the ret.time is within 0.15 min of each other
            # Checks the uv.count is at least 3 OR uv.mean < 1.5
            # If satisfied, the peak is declared a match and added to df
            
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
            
        }    else if (coculture$RetTime[cc.peak] < (con$RetTime[con.peak] + 0.15) && 
                      coculture$RetTime[cc.peak] > (con$RetTime[con.peak] -0.15) &&
                      final.count > 1 && abs(uv.mean) < 2) {  
            
            # Checks if the ret.time is within 0.15 min of each other
            # Checks the uv.count is at least 2 AND uv.mean < 2
            # Both uv.count and uv.mean are less strict requirement
            # But requiring satisfaction of both adds stringency
            # If satisfied, the peak is declared a match and added to df
            
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
            
        }    else if (coculture$RetTime[cc.peak] < (con$RetTime[con.peak] + 0.05) && 
                      coculture$RetTime[cc.peak] > (con$RetTime[con.peak] -0.05) &&
                      final.count > 1) {  
            
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
        }   else  { 
            
            # The peak is declared NOT to have a match in the control and given NA
            
            cc.df <- rbind(cc.df, 
                           c(PeakNo_CC = coculture$Peak[cc.peak], 
                             RetTime_CC = coculture$RetTime[cc.peak], 
                             PeakArea_CC = coculture$Area[cc.peak], 
                             PeakNo_con = NA, RetTime_con = NA, 
                             PeakArea_con = NA, UV_Count_con = NA, 
                             Subtracted_UV_Mean_con = NA, 
                             PeakRatio_con = NA))
        }   
        cc.peak <- cc.peak + 1
    }
    return(cc.df)
}

#############################################
# 6.ConConsolidator: Removes double peak matching to a coculture peak
#############################################

# The following code removes double assignments to a coculture peak
# from two different controls.

# An example is that perhaps peak 1 from the coculture matches
# to peak 1 in con1 and peak 2 in con2.

# This function will determine which control's matched peak matches
# closest and remove the assignment for the weakest match.

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
# 7.Metabolite Effect Characteriser: converting PeakRatio to a Factor
#############################################

# The following piece of code adds a column that categorises the peak 
# areas into suppressions, and enhancements

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
# MIMI_Part_1: The main working-function to compare the peak-matching and refinement
#############################################

MIMI <- function() {           
    
    #Will not work unless Interaction_Matrix object created by the user.
    
    matrix.total.rows <- nrow(Interaction_Matrix)
    matrix.row.no <- 1
    
    while (matrix.row.no <= matrix.total.rows) {
        
        con1.df <- as.character(Interaction_Matrix[matrix.row.no, 1])
        con2.df <- as.character(Interaction_Matrix[matrix.row.no, 2])
        cc.name <- as.character(Interaction_Matrix[matrix.row.no, 3])
        matrix.row.no <- matrix.row.no + 1
        print(cc.name)
        #Generates 3 dataframes using the ReadExcel function for the first
        #interction to be investigated from the Interaction_Matrix
        
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
        
        # Rewrites the tidied dataset to a csv file
        
        write.csv(refined.cc.df, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         cc.name, ".CSV"), row.names = FALSE)
    }
}

# Load in excel file named Interaction_Matrix and use MIMI() function:
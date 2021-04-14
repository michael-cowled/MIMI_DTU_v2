##########################################################
##########################################################
# Microbial Interaction Metabolite Integrator - Script 2 #
##########################################################
##########################################################

# COPYWRIGHT: © Macquarie University - Michael Cowled and contributors [2021]

# INPUT: The output file from Script_1 is used as the input file for Script_2.

# OUTPUT: A refined output file (.csv) is generated that remove doubley assigned 
# peaks, identifies non UVs and adds in peaks belonging to the control(s).


# The same R packages required of MIMI_Script_1 are required here.

#------------------------------------------------------------------------------#

# Functions to be pre-loaded prior to use of the main function, MIMI2()

# 1.SimpleEffectCategoriser - Uses the principles of the EffectCategoriser 
# function to categorise based on a single PeakRatio input.

# 2.IdentifyBadPeak: Decides on which of the doubley assigned peaks to remove

# 3.RemoveDoubleyAssignedPeaks -Checks for a peak in the control being assigned 
# to more than one peak in the coculture. Uses the number of matching UV maxima 
# and/or the subtracted UV spectra to make decisions as to which peak is a better 
# match.

# 4.PeakAssigner: Adds the assignment of the matched non UV peak from MatchNonUVs

# 5.MatchNonUVs - - Tentatively assigns matched peaks as non UVs 
# (or as distorted UVs) if matching the conditions.
 
# 6.InhibitionChecker - Checks if one culture is inhibited and deposits into a 
# separate output file.

# 7.RowBinder2: A variant of RowBinder used for adding missing control peaks

# 8.FindMissingcontrolPeaks - Adds in unassigned peaks from the controls to 
# provide a single, unified table.

#------------------------------------------------------------------------------#

# Note: Some functions as part of this script rely on functions from MIMI_Script_1.

# Reading in the functions:

#############################################
# 1.SimpleEffectCategoriser:
#############################################

# This is a simple version of the Effect categoriser with a single ratio input

SimpleEffectCategoriser <- function(ratio) {
    if (ratio > -100 && ratio <= 20) {
        Effect <- 2
    }   else if (ratio > -20 && ratio < 20) {
        Effect <- 3
    }   else if (ratio >= 20 && ratio < 100) {
        Effect <- 4
    }   else if (ratio >= 100) {
        Effect <- 5
    }
    return(Effect)
}

#############################################
# 2.IdentifyBadPeak: Decides on which of the doubley assigned peaks to remove
#############################################

IdentifyBadPeak <- function(double.peaks.subset) {
    
    if (double.peaks.subset[1, 10] > double.peaks.subset[2, 10]) {
        
        # Peaks are compared based on UV count first.
        # The peak with the lowest UV count is removed.
        
        bad.peak <- double.peaks.subset[2, 2]
    }   else if (double.peaks.subset[1, 10] < double.peaks.subset[2, 10]) {
        bad.peak <- double.peaks.subset[1, 2]
    }   else if (double.peaks.subset[1, 11] < double.peaks.subset[2, 11]) {
        
        # If the UV counts are equal the subtracted UV mean is compared.
        # The peak with the highest UV mean is removed.
        
        bad.peak <- double.peaks.subset[2, 2]
    }   else {
        bad.peak <- double.peaks.subset[1, 2]
    }
    return(bad.peak)
}


#############################################
# 3.RemoveDoubleyAssignedPeaks: Removes doubley-assigned peaks from a control.
#############################################

# An example is that perhaps the peaks 1 and 2 from the coculture match to
# twice to peak 1 in the control.

# This function will determine which peak from the control matches closest
# and remove the assignment for the weakest match.

RemoveDoubleyAssignedPeaks <- function(Interaction_Matrix) {
    
    logic.total.rows <- 1
    
    while (logic.total.rows > 0) {
        
        matrix.total.rows <- nrow(Interaction_Matrix)
        matrix.row.no <- 1
        logic.table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)),
                                c("cc.name"))
        
        # A df is set up to list the instances where double-peak matching occurs.
        
        while (matrix.row.no <= matrix.total.rows) {
            
            cc.name <- as.character(Interaction_Matrix[matrix.row.no,3])
            cc.name.list <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), 
                                     c("cc.name"))
            cc.name.list[1, 1] <- cc.name
            df.name <- 
                read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                                cc.name, ".csv"))
            df.name <- unite(df.name, Combined, c(Matched_con, PeakNo_con), 
                             sep = "-", remove = FALSE)
            df.name <- filter(df.name, Combined !="NA-NA")
            logic <- length(unique(df.name$Combined)) == nrow(df.name)
            
            # Compares the no. of unique peaks assigned in the control
            
            # If the no. of unique peaks differs to the no. of peaks in the
            # coculture, then it will be incorporated in the logic.table
            
            if (logic == TRUE) {
            } else {
                logic.table <- rbind(logic.table, cc.name.list)
            }
            matrix.row.no <- matrix.row.no +1
        }
        
        # Now to use the logic.table to read in the files that need fixing.
        
        # Note: Only 1 peak is fixed at a time, and so will go back through
        # and regenerate the logic.table and check if more peaks need fixing.
        
        logic.total.rows <- nrow(logic.table)
        logic.row.no <- 1
        
        while (logic.row.no <= logic.total.rows) {
            
            #Preprocessing code to read and manipulate the file of interest
            
            cc.name <- as.character(logic.table[logic.row.no, 1])
            double.peaks.df <- 
                read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                                cc.name, ".csv"))
            double.peaks.df <- unite(double.peaks.df, Combined, 
                                     c(Matched_con, PeakNo_con), sep = "-", 
                                     remove = FALSE)
            df.to.manipulate <- double.peaks.df
            double.peaks.df <- filter(double.peaks.df, Combined != "NA-NA")
            double.peaks.df$Duplicated <- duplicated(double.peaks.df$Combined)
            double.peaks.subset <- filter(double.peaks.df, Duplicated == TRUE)
            double.peaks.subset <- filter(double.peaks.df, 
                                          Combined == double.peaks.subset[1, 5])
            
            bad.peak <- IdentifyBadPeak(double.peaks.subset)
            
            df.to.manipulate[bad.peak, 5:12] <- NA
            double.peaks.df_removed <- select(df.to.manipulate, -Combined)
            write.csv(double.peaks.df_removed, 
                      paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                             cc.name, ".CSV"), row.names = FALSE)
            logic.row.no <- logic.row.no + 1
        }
    }
}

#############################################
# 4.PeakAssigner: Adds the assignment of the matched non UV peak from MatchNonUVs
#############################################

PeakAssigner <- function(df.name, cc.peak, con.name, con.peak, con, final.count,
                         ratio, Effect) {
    df.name$Matched_con[cc.peak] <-con.name
    df.name$PeakNo_con[cc.peak] <-con$Peak[con.peak]
    df.name$RetTime_con[cc.peak] <-con$RetTime[con.peak]
    df.name$PeakArea_con[cc.peak] <-con$Area[con.peak]
    df.name$UV_Count[cc.peak] <- final.count
    df.name$PeakRatio[cc.peak] <- ratio
    df.name$Metabolite_Effect[cc.peak] <- Effect
    return(df.name)
}

#############################################
# 5.MatchNonUVs: Further assigns peaks based on weaker criteria.
#############################################

MatchNonUVs <- function(Interaction_Matrix) {
    
    matrix.total.rows <- nrow(Interaction_Matrix)
    matrix.row.no <- 1
    
    while (matrix.row.no <= matrix.total.rows) {
        
        # Reads in the first coculture output file to be amended.
        # Reads in the the corresponding con files from raw NovaC.
        
        con1.name <- as.character(Interaction_Matrix[matrix.row.no, 1])
        con2.name <- as.character(Interaction_Matrix[matrix.row.no, 2])
        cc.name <- as.character(Interaction_Matrix[matrix.row.no, 3])
        df.name <- 
            read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                            cc.name, ".csv"))
        con1 <- as.data.frame(ReadExcel(con1.name))
        con2 <- as.data.frame(ReadExcel(con2.name))
        cc <- as.data.frame(ReadExcel(cc.name))
        
        matrix.row.no <- matrix.row.no + 1
        n <- nrow(df.name)
        cc.peak <- 1
        
        # 'cc.peak' corresponds to the peak no. to be matched in the coculture
        # This is a second set of peak matching that improves on the first round.
        
        while (cc.peak < n+1) {
            con1.peak <- which(abs(con1$RetTime-cc$RetTime[cc.peak]) ==
                                   min(abs(con1$RetTime-cc$RetTime[cc.peak])))
            con2.peak <- which(abs(con2$RetTime-cc$RetTime[cc.peak]) ==
                                   min(abs(con2$RetTime-cc$RetTime[cc.peak])))
            final.count.1 <- UVcheck(con1, cc, cc.peak,con1.peak)
            final.count.2 <- UVCheck(con2, cc, cc.peak,con2.peak)
            if (!is.na(df.name$Matched_con[cc.peak]) | 
                (any(df.name[, 7] ==con1$RetTime[con1.peak], na.rm = TRUE)) |
                (any(df.name[, 7] ==con2$RetTime[con2.peak], na.rm = TRUE))) {
                
                # Checks for peak matching already, and skips to the next peak.
                
            }   else if (cc$RetTime[cc.peak] < (con1$RetTime[con1.peak] + 0.02) && 
                         cc$RetTime[cc.peak] > (con1$RetTime[con1.peak] -0.02) &&
                         final.count.1 < 2)    {
                
                # Checks if the closest match incon1 satisfies this test.
                
                con.peak <-con1.peak
                ratio = (((cc[cc.peak, 3] -con1[con.peak, 3]) 
                          / con1[con.peak, 3]) * 100)
                Effect <- SimpleEffectCategoriser(ratio)
                df.name <- PeakAssigner(df.name, cc.peak, con1.name, con.peak, 
                                        con1, final.count.1, ratio, Effect)
                
            }   else if (cc$RetTime[cc.peak] < (con2$RetTime[con2.peak] + 0.02) && 
                         cc$RetTime[cc.peak] > (con2$RetTime[con2.peak] - 0.02) &&
                         final.count.2 < 2)    {
                con.peak <-con2.peak
                ratio = (((cc[cc.peak, 3] -con2[con.peak, 3]) / 
                              con2[con.peak, 3]) * 100)
                Effect <- SimpleEffectCategoriser(ratio)
                df.name <- PeakAssigner(df.name, cc.peak, con2.name, con.peak, 
                                        con2, final.count.2, ratio, Effect)
            }    
            cc.peak <- cc.peak + 1
        }
        write.csv(df.name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         cc.name, ".CSV"), row.names = FALSE)
    }
}

#############################################
# 6.InhibitionChecker: Verifies if one culture is inhibited
#############################################

InhibitionChecker <- function(df.name, inhibition.df, cc.name) {
    
    #From opened output file
    
    match <- as.data.frame(table(df.name$Matched_con))
    names(match)[2] <- 'match'
    unmatch <- as.data.frame(table(df.name$Sample_Ref))
    names(unmatch)[2] <- 'unmatch'
    combined <- merge(match, unmatch) %>%
        mutate(ratio = unmatch / match)
    
    # Then verify whether inhibition of control indicated
    
    if (nrow(combined) == 0) {
        combined <- match
        combined$unmatch[[1]] <- 0
        combined$ratio[[1]] <- 0
        combined$Inhibition[[1]] <- TRUE
    } else if (nrow(combined) == 1) {
        combined$Inhibition[[1]] <- TRUE
    }   else if (combined$match[1] <= 1 && combined$ratio[1] >= 3) {
        combined$Inhibition[[1]] <- TRUE
        combined$Inhibition[[2]] <- FALSE
    }   else if (combined$match[2] <= 1 && combined$ratio[2] >= 3) {
        combined$Inhibition[[2]] <- TRUE
        combined$Inhibition[[1]] <- FALSE
    }   else {
        combined$Inhibition[[2]] <- FALSE
        combined$Inhibition[[1]] <- FALSE
    }
    
    combined <- filter(combined, Inhibition == TRUE)
    
    if (nrow(combined) > 0) {
        temp <- setNames(data.frame(matrix(ncol = 3, nrow = 1)),
                         c("cc.name", "Inhibition", "Dominating_Culture"))
        combined$cc.name <- cc.name
        combined$Var1 <- as.character(combined$Var1)
        temp$cc.name <- combined[1, 6]
        temp$Inhibition <- combined[1, 5]
        temp$Dominating_Culture <- combined[1, 1]
        inhibition.df <- rbind(inhibition.df, temp)
    }
    return(inhibition.df)
}

#############################################
# 7.RowBinder2: A variant of RowBinder used for adding missing control peaks
#############################################

RowBinder2 <- function(df.name, con.name, con, cc.peak) {
    
    df.name <- rbind(df.name, 
                     c(Sample_Ref = con.name, PeakNo_CC = NA, 
                       RetTime_CC = NA, PeakArea_CC = NA, 
                       Combined = NA, Matched_con = NA, 
                       PeakNo_con =con$Peak[cc.peak], 
                       RetTime_con =con$RetTime[cc.peak], 
                       PeakArea_con =con$Area[cc.peak], 
                       UV_Count = NA, Subtracted_UV_Mean = NA, 
                       PeakRatio = -100, Metabolite_Effect = 1))
    return(df.name)
}

#############################################
# 8.FindMissingcontrolPeaks: Adds in the unassigned peaks from the control(s)
#############################################

FindMissingcontrolPeaks <- function(Interaction_Matrix) {
    
    # Makes a new table that is used in InhibitionChecker function:
    inhibition.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                              c("cc.name", "Inhibition", "Inhibited_Culture"))
    
    matrix.total.rows <- nrow(Interaction_Matrix)
    matrix.row.no <- 1
    
    while (matrix.row.no <= matrix.total.rows) {
        
        # Reads in the first coculture output file to be amended.
        # Reads in the the corresponding con files from raw NovaC.
        
        con1.name <- as.character(Interaction_Matrix[matrix.row.no, 1])
        con2.name <- as.character(Interaction_Matrix[matrix.row.no, 2])
        cc.name <- as.character(Interaction_Matrix[matrix.row.no, 3])
        df.name <- 
            read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                            cc.name, ".csv"))
        
        df.name$Sample_Ref <- as.character(df.name$Sample_Ref)
        df.name <- unite(df.name, Combined, c(Matched_con, PeakNo_con), 
                         sep = "-", remove = FALSE)
        con1 <- as.data.frame(ReadExcel(con1.name))
        con2 <- as.data.frame(ReadExcel(con2.name))
        
        matrix.row.no <- matrix.row.no + 1
        n <- nrow(con1)
        cc.peak <- 1
        
        # Sequentially checkscon1 for peak 'cc.peak' in coculture output file
        
        while (cc.peak <= n) {
            if (any(df.name[, 5] == paste0(con1.name, "-", cc.peak), na.rm = TRUE)) {
            }   else if (any(df.name[, 5] != paste0(con1.name, "-", cc.peak), 
                             na.rm = TRUE)) {
                df.name <- RowBinder2(df.name, con1.name, con1, cc.peak)
            }
            cc.peak <- cc.peak + 1 
        }
        
        n <- nrow(con2)
        cc.peak <- 1
        
        while (cc.peak <= n) {
            if (any(df.name[, 5] == paste0(con2.name, "-", cc.peak), na.rm = TRUE)) {
            }   else if (any(df.name[, 5] != paste0(con2.name, "-", cc.peak), 
                             na.rm = TRUE)) {
                df.name <- RowBinder2(df.name, con2.name, con2, cc.peak)
            }   
            cc.peak <- cc.peak + 1
        }
        df.name <- select(df.name, -Combined)
        
        # Space to include function to identify inhibition
        inhibition.df <- InhibitionChecker(df.name, inhibition.df, 
                                           cc.name)
        write.csv(df.name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         cc.name, ".CSV"), row.names = FALSE)
    }
    inhibition.df <- transform(inhibition.df, Inhibition = as.logical(Inhibition))
    write.csv(inhibition.df, 
              paste0("Testing Broad-Scale Interactions/OutputFiles/inhibition.df.CSV"), 
              row.names = FALSE)
}

#############################################
# MIMI_Part_2: Carries out the main dereplication processes
#############################################

# Carries out the last two functions:
# 1.RemoveDoubleyAssignedPeaks: Multiple peaks in a coculture matched to the same
# unique peak of a control
# 2.FindMissingcontrolPeaks: Unique peaks from control(s) not matched to a peak
# in the coculture, are added into a single, unified df

MIMI2 <- function() {  
    
    print("Initiating RemoveDoubleyAssignedPeaks")
    RemoveDoubleyAssignedPeaks(Interaction_Matrix)
    print("Initiating MatchNonUVs")
    MatchNonUVs(Interaction_Matrix)
    print("Initiating FindMissingcontrolPeaks")
    FindMissingcontrolPeaks(Interaction_Matrix)
    print("MIMI2 completed.")
    
}
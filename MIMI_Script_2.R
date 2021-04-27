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

# 2.DoublePeakChecker: Finds and reports instances of double peak matching

# 3.IdentifyBadPeak: Decides on which of the doubley assigned peaks to remove

# 4.RemoveDoubleyAssignedPeaks -Checks for a peak in the control being assigned 
# to more than one peak in the coculture. Uses the number of matching UV maxima 
# and/or the subtracted UV spectra to make decisions as to which peak is a better 
# match.

# 5.PeakAssigner: Adds the assignment of the matched non UV peak from MatchNonUVs

# 6.MatchNonUVs - - Tentatively assigns matched peaks as non UVs 
# (or as distorted UVs) if matching the conditions.

# 7.RowBinder2: A variant of RowBinder used for adding missing control peaks

# 8.FindMissingcontrolPeaks - Adds in unassigned peaks from the controls to 
# provide a single, unified table.

#------------------------------------------------------------------------------#

# Note: Some functions as part of this script rely on functions from MIMI_Script_1.

# Reading in the functions:

#############################################
# 1.SimpleEffectCategoriser
#############################################

# Characterises the ratios of peak area changes to a categorical factor.
# Only returns the "Effect" variable, rather than updating the df.

# Args:
# ratio is calculated from peak area previously

# Returns:
# Effect - a categorised effect on the metabolite based on ratio range.

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
# 2.DoublePeakChecker
#############################################

# Finds and reports instances of double peak matching in a "peak.misassignment.list"

# Args:
# Interaction_Matrix to determine output files to check

# Returns:
# peak.misassignment.list - lists the output files containing doubley assigned peaks.

DoublePeakChecker <- function(Interaction_Matrix) {
    
    peak.misassignment.list <- setNames(data.frame(matrix(ncol = 1, nrow = 0)),
                            c("cc.name"))
    
    matrix.total.rows <- nrow(Interaction_Matrix)
    matrix.row.no <- 1
    
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
        # coculture, then it will be incorporated in the peak.misassignment.list
        
        if (logic == TRUE) {
        } else {
            peak.misassignment.list<- rbind(peak.misassignment.list, cc.name.list)
        }
        matrix.row.no <- matrix.row.no +1
    }
    return(peak.misassignment.list)
}

#############################################
# 3.IdentifyBadPeak
#############################################

# Decides on which of the doubley assigned peaks to remove

# Args:
# double.peaks.subset is the row subset of two peaks in the coculture which have 
# been assigned to more than one peak in the control.

# Returns:
# bad.peak - the peak to be removed.

IdentifyBadPeak <- function(double.peaks.subset) {
    
    if (double.peaks.subset[1, 11] > double.peaks.subset[2, 11]) {
        
        # Peaks are compared based on UV count first.
        # The peak with the lowest UV count is removed.
        
        bad.peak <- double.peaks.subset[2, 2]
    }   else if (double.peaks.subset[1, 11] < double.peaks.subset[2, 11]) {
        bad.peak <- double.peaks.subset[1, 2]
    }   else if (double.peaks.subset[1, 12] < double.peaks.subset[2, 12]) {
        
        # If the UV counts are equal the subtracted UV mean is then compared.
        # The peak with the highest UV mean is removed.
        
        bad.peak <- double.peaks.subset[2, 2]
    }   else {
        bad.peak <- double.peaks.subset[1, 2]
    }
    return(bad.peak)
}

#############################################
# 4.RemoveDoubleyAssignedPeaks
#############################################

# Removes doubley-assigned peaks from a control.

# An example is that perhaps the peaks 1 and 2 from the coculture match to
# twice to peak 1 in the control.

# This function will determine which peak from the control matches closest
# and remove the assignment for the weakest match.

# Args:
# Interaction_Matrix to be fed into the DoublePeakChecker function.

# Returns:
# double.peaks.df_removed - an updated df with one instance of double peak
# matching fixed; but is reiterated through function.

RemoveDoubleyAssignedPeaks <- function(Interaction_Matrix) {
    
    logic.total.rows <- 1
    
    while (logic.total.rows > 0) {
        
        peak.misassignment.list <- DoublePeakChecker(Interaction_Matrix)
        
        # Now to use the peak.misassignment.list to read in the files that need fixing.
        
        # Note: Only 1 peak is fixed at a time, and so will go back through
        # and regenerate the peak.misassignment.list and check if more peaks need fixing.
        
        if (is.null(nrow(peak.misassignment.list))) {
            logic.total.rows <- 0
        }   else {
            logic.total.rows <- nrow(peak.misassignment.list)
        }
    }
    logic.row.no <- 1
    
    while (logic.row.no <= logic.total.rows) {
        
        #Preprocessing code to read and manipulate the file of interest
        
        cc.name <- as.character(peak.misassignment.list[logic.row.no, 1])
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
                                      Combined == double.peaks.subset[1, 6])
        
        bad.peak <- IdentifyBadPeak(double.peaks.subset)
        
        df.to.manipulate[bad.peak, 5:13] <- NA
        double.peaks.df_removed <- select(df.to.manipulate, -Combined)
        write.csv(double.peaks.df_removed, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         cc.name, ".CSV"), row.names = FALSE)
        logic.row.no <- logic.row.no + 1
        return(double.peaks.df_removed)
    }
}

#############################################
# 5.PeakAssigner
#############################################

# Adds the assignment of the matched non UV peak from MatchNonUVs

# Args:
# df.name is the output file to be edited.
# cc.peak is the peak number of the coculture to be matched.
# con.name is the name of the control to which the matched peak corresponds.
# con.peak is the peak number of the control to with the matched peak corresponds.
# con is the df generated for the control to which the matched peak corresponds.
# final.count is the number of matching uv maxima between cc and con peaks.
# ratio is the peak area ratio for the matching cc and con peaks.
# Effect is the categorised effect on the metabolite as stipulated by the ratio.

# Returns:
# df.name - the updated df with the new peak assigned.

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
# 6.MatchNonUVs
#############################################

# Further assigns peaks based on weaker criteria.

# Args: 
# cc.peak is the peak number of the coculture to be matched.
# con1 is the df generated for con1.
# con2 is the df generated for con2.
# cc is the df generated for cc.
# df.name is the df to be manipulated.
# con1.name is the name of con1.
# con2.name is the name of con2.

# Returns:
# df.name - the updated df with the non-UV peak assigned.

MatchNonUVs <- function(cc.peak, con1, con2, cc, df.name, con1.name, con2.name) {
    
    n <- nrow(df.name)
    
    while (cc.peak < n + 1) {
        con1.peak <- which(abs(con1$RetTime-cc$RetTime[cc.peak]) ==
                               min(abs(con1$RetTime-cc$RetTime[cc.peak])))
        con2.peak <- which(abs(con2$RetTime-cc$RetTime[cc.peak]) ==
                               min(abs(con2$RetTime-cc$RetTime[cc.peak])))
        final.count.1 <- CheckUVCount(con1, cc, cc.peak,con1.peak)
        final.count.2 <- CheckUVCount(con2, cc, cc.peak,con2.peak)
        if (!is.na(df.name$Matched_con[cc.peak]) | 
            (any(df.name[, 7] ==con1$RetTime[con1.peak], na.rm = TRUE)) |
            (any(df.name[, 7] ==con2$RetTime[con2.peak], na.rm = TRUE))) {
            
            # Checks for peak matching already, and skips to the next peak.
            
        }   else if (cc$RetTime[cc.peak] < (con1$RetTime[con1.peak] + 0.02) && 
                     cc$RetTime[cc.peak] > (con1$RetTime[con1.peak] -0.02) &&
                     final.count.1 < 2)    {
            
            # Checks if the closest match in con1 satisfies this test.
            
            ratio <- CalcRatio(cc, cc.peak, con1, con1.peak)
            Effect <- SimpleEffectCategoriser(ratio)
            df.name <- PeakAssigner(df.name, cc.peak, con1.name, con1.peak, 
                                    con1, final.count.1, ratio, Effect)
            
        }   else if (cc$RetTime[cc.peak] < (con2$RetTime[con2.peak] + 0.02) && 
                     cc$RetTime[cc.peak] > (con2$RetTime[con2.peak] - 0.02) &&
                     final.count.2 < 2)    {
            ratio <- CalcRatio(cc, cc.peak, con2, con2.peak)
            Effect <- SimpleEffectCategoriser(ratio)
            df.name <- PeakAssigner(df.name, cc.peak, con2.name, con2.peak, 
                                    con2, final.count.2, ratio, Effect)
        }    
        cc.peak <- cc.peak + 1
    }
    return(df.name)
}

#############################################
# 7.RowBinder2
#############################################

# A variant of the RowBinder function; to be used for adding missing control peaks
# Notably there is NA values associated with the cc as no match found.

# Args:
# df.name is the df to be manipulated.
# con.name is the name of the control to which the matched peak corresponds.
# con is the df generated for the control to which the matched peak corresponds. 
# cc.peak is the peak number of the coculture to be matched.

# Returns:
# df.name - the updated df with the new peak assigned.

RowBinder2 <- function(df.name, con.name, con, cc.peak) {
    
    df.name <- rbind(df.name, 
                     c(Sample_Ref = con.name, PeakNo_CC = NA, 
                       RetTime_CC = NA, PeakArea_CC = NA, PercArea = NA, 
                       Combined = NA, Matched_con = NA, 
                       PeakNo_con =con$Peak[cc.peak], 
                       RetTime_con =con$RetTime[cc.peak], 
                       PeakArea_con =con$Area[cc.peak], 
                       UV_Count = NA, Subtracted_UV_Mean = NA, 
                       PeakRatio = -100, Metabolite_Effect = 1))
    return(df.name)
}

#############################################
# 8.FindMissingcontrolPeaks
#############################################

# Adds in the unassigned peaks from the control(s).

# Args:
# df.name is the df to be manipulated.
# con1.name is the name of con1.
# con1 is the df generated for con1.
# con2.name is the name of con2. 
# con2 is the df generated for con2.

#Returns:
# # df.name - the updated df with the missing control peaks assigned.

FindMissingcontrolPeaks <- function(df.name, con1.name, con1, con2.name, con2) {
    
    n <- nrow(con1)
    cc.peak <- 1
    
    # Sequentially checks con1 for peak 'cc.peak' in coculture output file
    
    while (cc.peak <= n) {
        if (any(df.name[, 6] == paste0(con1.name, "-", cc.peak), na.rm = TRUE)) {
        }   else if (any(df.name[, 6] != paste0(con1.name, "-", cc.peak), 
                         na.rm = TRUE)) {
            df.name <- RowBinder2(df.name, con1.name, con1, cc.peak)
        }
        cc.peak <- cc.peak + 1 
    }
    
    n <- nrow(con2)
    cc.peak <- 1
    
    while (cc.peak <= n) {
        if (any(df.name[, 6] == paste0(con2.name, "-", cc.peak), na.rm = TRUE)) {
        }   else if (any(df.name[, 6] != paste0(con2.name, "-", cc.peak), 
                         na.rm = TRUE)) {
            df.name <- RowBinder2(df.name, con2.name, con2, cc.peak)
        }   
        cc.peak <- cc.peak + 1
    }
    df.name <- select(df.name, -Combined)
    return(df.name)
}

#############################################
# MIMI_Part_2: Carries out the main dereplication processes
#############################################

# Microbial Interaction Metabolite Integrator.
# The secondary function to dereplicate peak matching, add non-UVs and missing
# control peaks.

# Args:
# Interaction_Matrix is a user defined matrix contraining the column names:
# CON1
# CON2	
# Coculture
# Note: Output files generated from MIMI script 1 should be present for each cc.
# Also note: NovaC files (a refined version of the raw HPLC data) should be 
# present for each control and coculture.

# Returns:
# Finalised output files derived from df.name which summarise all matched peaks
# in each coculture and its corresponding controls.

MIMI2 <- function() {  
    
    print("Initiating RemoveDoubleyAssignedPeaks")
    RemoveDoubleyAssignedPeaks(Interaction_Matrix)
    
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
        cc.peak <- 1
        
        df.name <- MatchNonUVs(cc.peak, con1, con2, cc, df.name, con1.name, con2.name)
        df.name$Sample_Ref <- as.character(df.name$Sample_Ref)
        df.name <- unite(df.name, Combined, c(Matched_con, PeakNo_con), 
                         sep = "-", remove = FALSE)
        df.name <- FindMissingcontrolPeaks(df.name, con1.name, con1, con2.name, con2)
        
        write.csv(df.name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         cc.name, ".CSV"), row.names = FALSE)
    }
    print("MIMI2 completed.")
}
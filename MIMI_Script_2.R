#############################################
#############################################
#Microbial Interaction Metabolite Integrator#
#############################################
#############################################

#A Matrix defined as Interaction_Matrix is required to be made with the THREE
#samples to be compared.

#R packages required to be loaded in:

library(dplyr)
library(tidyr)
library(readxl)

#Functions to be pre-loaded prior to use of the main function, MIMI()
#1.SimpleEffectCategoriser
#2.DoublePeakRemover
#3.NonUVPeakMatcher
#4.InhibitionChecker
#5.MissingControlPeaks

#Note: Some functions as part of this script rely on functions from MIMI_Script_1.

#Reading in the functions:

#############################################
#1.Simple Metabolite Effect Characteriser:
#############################################

#This is a simple version of the effect categoriser with a single ratio input

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
#2.DoublePeakRemover: Removes doubley-assigned peaks from a control.
#############################################

#An example is that perhaps the peaks 1 and 2 from the coculture match to
#twice to peak 1 in the control.

#This function will determine which peak from the control matches closest
#and remove the assignment for the weakest match.

DoublePeakRemover <- function(Interaction_Matrix) {
    
    Logic.Total.Rows <- 1
    
    while (Logic.Total.Rows > 0) {
        
        Matrix.Total.Rows <- nrow(Interaction_Matrix)
        Matrix.Row.no <- 1
        Logic.Table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)),
                                c("CC.Name"))
        
        #A df is set up to list the instances where double-peak matching occurs.
        
        while (Matrix.Row.no <= Matrix.Total.Rows) {
            
            CC.Name <- as.character(Interaction_Matrix[Matrix.Row.no,3])
            CC.Name.List <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), 
                                        c("CC.Name"))
            CC.Name.List[1,1] <- CC.Name
            df.Name <- 
                read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                                CC.Name, ".csv"))
            df.Name <- unite(df.Name, Combined, c(Matched_CON, PeakNo_CON), 
                             sep = "-", remove = FALSE)
            df.Name <- filter(df.Name, Combined !="NA-NA")
            logic <- length(unique(df.Name$Combined)) == nrow(df.Name)
            
            #Compares the no. of unique peaks assigned in the control
            
            #If the no. of unique peaks differs to the no. of peaks in the
            #coculture, then it will be incorporated in the Logic.Table
            
            if (logic == TRUE) {
            } else {
                Logic.Table <- rbind(Logic.Table, CC.Name.List)
            }
            Matrix.Row.no <- Matrix.Row.no +1
        }
        
        #Now to use the Logic.Table to read in the files that need fixing.
        
        #Note: Only 1 peak is fixed at a time, and so will go back through
        #and regenerate the Logic.Table and check if more peaks need fixing.
        
        Logic.Total.Rows <- nrow(Logic.Table)
        Logic.Row.no <- 1

        while (Logic.Row.no <= Logic.Total.Rows) {
            
            #Preprocessing code to read and manipulate the file of interest
            
            CC.Name <- as.character(Logic.Table[Logic.Row.no,1])
            df.with.double.peaks <- 
                read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                                CC.Name, ".csv"))
            df.with.double.peaks <- unite(df.with.double.peaks, Combined, c(Matched_CON, PeakNo_CON), 
                             sep = "-", remove = FALSE)
            df.to.manipulate <- df.with.double.peaks
            df.with.double.peaks <- filter(df.with.double.peaks, Combined != "NA-NA")
            df.with.double.peaks$Duplicated <- duplicated(df.with.double.peaks$Combined)
            df.to.compare.double.peaks <- filter(df.with.double.peaks, Duplicated == TRUE)
            df.to.compare.double.peaks <- filter(df.with.double.peaks, Combined == df.to.compare.double.peaks[1, 5])
            
            if (df.to.compare.double.peaks[1,10] > df.to.compare.double.peaks[2,10]) {
                
                #Peaks are compared based on UV count first.
                #The peak with the lowest UV count is removed.
                Bad.Peak <- df.to.compare.double.peaks[2,2]
            }   else if (df.to.compare.double.peaks[1,10] < df.to.compare.double.peaks[2,10]) {
                Bad.Peak <- df.to.compare.double.peaks[1,2]
            }   else if (df.to.compare.double.peaks[1,11] < df.to.compare.double.peaks[2,11]) {
                
                #If the UV counts are equal the subtracted UV mean is compared.
                #The peak with the highest UV mean is removed.
                
                Bad.Peak <- df.to.compare.double.peaks[2,2]
            }   else {
                Bad.Peak <- df.to.compare.double.peaks[1,2]
            }
            df.to.manipulate[Bad.Peak, 5:12] <- NA
            df.with.double.peaks_removed <- select(df.to.manipulate, -Combined)
            write.csv(df.with.double.peaks_removed, 
                      paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                             CC.Name, ".CSV"), row.names = FALSE)
            Logic.Row.no <- Logic.Row.no + 1
        }
    }
}

#############################################
#3.NON_UV Peak Matcher: Further assigns peaks based on weaker criteria.
#############################################

NonUVMatcher <- function(Interaction_Matrix) {
    
    Matrix.Total.Rows <- nrow(Interaction_Matrix)
    Matrix.Row.no <- 1
    
    while (Matrix.Row.no <= Matrix.Total.Rows) {
        
        #Reads in the first coculture output file to be amended.
        #Reads in the the corresponding CON files from raw NovaC.
        
        CON1.Name <- as.character(Interaction_Matrix[Matrix.Row.no,1])
        CON2.Name <- as.character(Interaction_Matrix[Matrix.Row.no,2])
        CC.Name <- as.character(Interaction_Matrix[Matrix.Row.no,3])
        df.Name <- 
            read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                            CC.Name, ".csv"))
        CON1 <- as.data.frame(ReadExcel(CON1.Name))
        CON2 <- as.data.frame(ReadExcel(CON2.Name))
        Coculture <- as.data.frame(ReadExcel(CC.Name))
        
        Matrix.Row.no <- Matrix.Row.no + 1
        n <- nrow(df.Name)
        CC.peak <- 1
        
        #'CC.peak' corresponds to the peak no. to be matched in the coculture
        #This is a second set of peak matching that improves on the first round.
        
        while (CC.peak < n+1) {
            CON1.peak <- which(abs(CON1$RetTime-Coculture$RetTime[CC.peak]) ==
                            min(abs(CON1$RetTime-Coculture$RetTime[CC.peak])))
            CON2.peak <- which(abs(CON2$RetTime-Coculture$RetTime[CC.peak]) ==
                            min(abs(CON2$RetTime-Coculture$RetTime[CC.peak])))
            Final.Count.1 <- UVcheck(CON1, Coculture, CC.peak, CON1.peak)
            Final.Count.2 <- UVCheck(CON2, Coculture, CC.peak, CON2.peak)
            if (!is.na(df.Name$Matched_CON[CC.peak]) | 
                (any(df.Name[,7] == CON1$RetTime[CON1.peak], na.rm = TRUE)) |
                (any(df.Name[,7] == CON2$RetTime[CON2.peak], na.rm = TRUE))) {
                
                #Checks for peak matching already, and skips to the next peak.

            }   else if (Coculture$RetTime[CC.peak] < (CON1$RetTime[CON1.peak] + 0.05) && 
                         Coculture$RetTime[CC.peak] > (CON1$RetTime[CON1.peak] -0.05) &&
                         Final.Count.1 < 2)    {
                
                #Checks if the closest match in CON1 satisfies this test.
                
                CON.peak <- CON1.peak
                ratio = (((Coculture[CC.peak,3] - CON1[CON.peak,3])/CON1[CON.peak,3])*100)
                Effect <- SimpleEffectCategoriser(ratio)
                
                #Assignments of the matched peak
                
                df.Name$Matched_CON[CC.peak] <- CON1.Name
                df.Name$PeakNo_CON[CC.peak] <- CON1$Peak[CON.peak]
                df.Name$RetTime_CON[CC.peak] <- CON1$RetTime[CON.peak]
                df.Name$PeakArea_CON[CC.peak] <- CON1$Area[CON.peak]
                df.Name$UV_Count[CC.peak] <- Final.Count.1
                df.Name$PeakRatio[CC.peak] <- ratio
                df.Name$Metabolite_Effect[CC.peak] <- Effect
                
            }   else if (Coculture$RetTime[CC.peak] < (CON2$RetTime[CON2.peak] + 0.05) && 
                         Coculture$RetTime[CC.peak] > (CON2$RetTime[CON2.peak] -0.05) &&
                         Final.Count.2 < 2)    {
                CON.peak <- CON2.peak
                ratio = (((Coculture[CC.peak,3] - CON2[CON.peak,3])/CON2[CON.peak,3])*100)
                Effect <- SimpleEffectCategoriser(ratio)
                df.Name$Matched_CON[CC.peak] <- CON2.Name
                df.Name$PeakNo_CON[CC.peak] <- CON2$Peak[CON.peak]
                df.Name$RetTime_CON[CC.peak] <- CON2$RetTime[CON.peak]
                df.Name$PeakArea_CON[CC.peak] <- CON2$Area[CON.peak]
                df.Name$UV_Count[CC.peak] <- Final.Count.2
                df.Name$PeakRatio[CC.peak] <- ratio
                df.Name$Metabolite_Effect[CC.peak] <- Effect
            }    
            CC.peak <- CC.peak + 1
        }
        write.csv(df.Name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         CC.Name, ".CSV"), row.names = FALSE)
    }
}

#############################################
#4.Inhibition Checker: Verifies if one culture is inhibited
#############################################

InhibitionChecker <- function(df.Name, Inhibition.df, CC.Name) {
    
    #From an opened output file (might be possible just downstream of missing control peaks)
    match <- as.data.frame(table(df.Name$Matched_CON))
    names(match)[2] <- 'match'
    unmatch <- as.data.frame(table(df.Name$Sample_Ref))
    names(unmatch)[2] <- 'unmatch'
    combined <- merge(match, unmatch) %>%
        mutate(ratio = unmatch / match)
    
    #Then verify whether inhibition of control indicated
    
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
                         c("CC.Name", "Inhibition", "Dominating_Culture"))
        combined$CC.Name <- CC.Name
        combined$Var1 <- as.character(combined$Var1)
        temp$CC.Name <- combined[1,6]
        temp$Inhibition <- combined[1,5]
        temp$Dominating_Culture <- combined[1,1]
        Inhibition.df <- rbind(Inhibition.df, temp)
    }
    return(Inhibition.df)
}

#############################################
#5.Missing Control Peaks: Adds in the unassigned peaks from the control(s)
#############################################

MissingControlPeaks <- function(Interaction_Matrix) {
    
    #Makes a new table that is used in InhibitionChecker function:
    Inhibition.df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                              c("CC.Name", "Inhibition", "Inhibited_Culture"))
    
    Matrix.Total.Rows <- nrow(Interaction_Matrix)
    Matrix.Row.no <- 1
    
    while (Matrix.Row.no <= Matrix.Total.Rows) {
        
        #Reads in the first coculture output file to be amended.
        #Reads in the the corresponding CON files from raw NovaC.
        
        CON1.Name <- as.character(Interaction_Matrix[Matrix.Row.no,1])
        CON2.Name <- as.character(Interaction_Matrix[Matrix.Row.no,2])
        CC.Name <- as.character(Interaction_Matrix[Matrix.Row.no,3])
        df.Name <- 
            read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                            CC.Name, ".csv"))
        
        df.Name$Sample_Ref <- as.character(df.Name$Sample_Ref)
        df.Name <- unite(df.Name, Combined, c(Matched_CON, PeakNo_CON), 
                         sep = "-", remove = FALSE)
        CON1 <- as.data.frame(ReadExcel(CON1.Name))
        CON2 <- as.data.frame(ReadExcel(CON2.Name))
        
        Matrix.Row.no <- Matrix.Row.no + 1
        n <- nrow(CON1)
        CC.peak <- 1
        
        #Sequentially checks CON1 for peak 'CC.peak' in coculture output file

        while (CC.peak <= n) {
            if (any(df.Name[,5] == paste0(CON1.Name, "-", CC.peak), na.rm = TRUE)) {
            }   else if (any(df.Name[,5] != paste0(CON1.Name, "-", CC.peak), 
                             na.rm = TRUE)) {
                df.Name <- rbind(df.Name, 
                                 c(Sample_Ref = CON1.Name, PeakNo_CC = NA, 
                                   RetTime_CC = NA, PeakArea_CC = NA, 
                                   Combined = NA, MAtched_CON = NA, 
                                   PeakNo_CON = CON1$Peak[CC.peak], 
                                   RetTime_CON = CON1$RetTime[CC.peak], 
                                   PeakArea_CON = CON1$Area[CC.peak], 
                                   UV_Count = NA, Subtracted_UV_Mean = NA, 
                                   PeakRatio = -100, Metabolite_Effect = 1))
            }
            CC.peak <- CC.peak + 1 
        }
        
        n <- nrow(CON2)
        CC.peak <- 1
        
        while (CC.peak <= n) {
            if (any(df.Name[,5] == paste0(CON2.Name, "-", CC.peak), na.rm = TRUE)) {
            }   else if (any(df.Name[,5] != paste0(CON2.Name, "-", CC.peak), 
                             na.rm = TRUE)) {
                df.Name <- rbind(df.Name, 
                                 c(Sample_Ref = CON2.Name, PeakNo_CC = NA, 
                                   RetTime_CC = NA, PeakArea_CC = NA,
                                   Combined = NA, MAtched_CON = NA, 
                                   PeakNo_CON = CON2$Peak[CC.peak], 
                                   RetTime_CON = CON2$RetTime[CC.peak], 
                                   PeakArea_CON = CON2$Area[CC.peak], 
                                   UV_Count = NA, Subtracted_UV_Mean = NA, 
                                   PeakRatio = -100, Metabolite_Effect = 1))
            }   
            CC.peak <- CC.peak + 1
        }
        df.Name <- select(df.Name, -Combined)
        
        #Space to include function to identify inhibition
        Inhibition.df <- InhibitionChecker(df.Name, Inhibition.df, 
                                            CC.Name)
        write.csv(df.Name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         CC.Name, ".CSV"), row.names = FALSE)
    }
    Inhibition.df <- transform(Inhibition.df, Inhibition = as.logical(Inhibition))
    write.csv(Inhibition.df, 
              paste0("Testing Broad-Scale Interactions/OutputFiles/Inhibition.df.CSV"), 
              row.names = FALSE)
}

#############################################
#MIMI_Part_2: Carries out the main dereplication processes
#############################################

#Carries out the last two functions:
#1.DoublePeakRemover: Multiple peaks in a coculture matched to the same
#unique peak of a control
#2.MissingControlPeaks: Unique peaks from control(s) not matched to a peak
#in the coculture, are added into a single, unified df

MIMI2 <- function() {  
    
print("Initiating DoublePeakRemover")
DoublePeakRemover(Interaction_Matrix)
print("Initiating NonUVMatcher")
NonUVMatcher(Interaction_Matrix)
print("Initiating MissingControlPeaks")
MissingControlPeaks(Interaction_Matrix)
print("MIMI2 completed.")

}
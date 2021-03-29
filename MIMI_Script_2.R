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
#1.Simple_Effect_Categoriser
#1.Double_Peak_Remover
#2.NON_UV_Peak_Matcher
#3.Inhibition_Checker
#4.Missing_Control_Peaks

#Reading in the functions:

#############################################
#1.Simple Metabolite Effect Characteriser:
#############################################

#This is a simple version of the effect categoriser with a single ratio input

Simple_Effect_Categoriser <- function(ratio) {
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
#2.Double Peak Remover: Removes doubley-assigned peaks from a control.
#############################################

#An example is that perhaps the peaks 1 and 2 from the coculture match
#twice to peak 1 in the control.

#This function will determine which peak from the control matches closest
#and remove the assignment for the weakest match.

Double_Peak_Remover <- function(Interaction_Matrix) {
    
    #Arbitrarily sets a positive value to the variable.
    Logic_TotalRows <- 1
    
    while (Logic_TotalRows > 0) {
        
        Matrix_TotalRows <- nrow(Interaction_Matrix)
        Matrix_Row_No <- 1
        Logic_Table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)),
                                c("Coculture_Name"))
        
        #A df is set up to list the instances where double-peak matching occurs.
        
        while (Matrix_Row_No <= Matrix_TotalRows) {
            
            Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
            Coculture_Name2 <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), 
                                        c("Coculture_Name"))
            Coculture_Name2[1,1] <- Coculture_Name
            df_Name <- 
                read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                                Coculture_Name, ".csv"))
            df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), 
                             sep = "-", remove = FALSE)
            df_Name <- filter(df_Name, Combined !="NA-NA")
            logic <- length(unique(df_Name$Combined)) == nrow(df_Name)
            
            #Compares the no. of unique peaks assigned in the control
            
            #If the no. of unique peaks differs to the no. of peaks in the
            #coculture, then it will be incorporated in the logic_table
            
            if (logic == TRUE) {
                Matrix_Row_No <- Matrix_Row_No +1
            } else {
                Logic_Table <- rbind(Logic_Table, Coculture_Name2)
                Matrix_Row_No <- Matrix_Row_No +1
            }
        }
        
        #Now to use the logic_table to read in the files that need fixing.
        
        #Note: Only 1 peak is fixed at a time, and so will go back through
        #and regenerate the logic_table and check if more peaks need fixing.
        
        Logic_TotalRows <- nrow(Logic_Table)
        Logic_Row_No <- 1
        
        while (Logic_Row_No <= Logic_TotalRows) {
            
            #Preprocessing code to read and manipulate the file of interest
            
            Coculture_Name <- as.character(Logic_Table[Logic_Row_No,1])
            df_Name <- 
                read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                                Coculture_Name, ".csv"))
            df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), 
                             sep = "-", remove = FALSE)
            df_Name2 <- filter(df_Name, Combined != "NA-NA")
            df_Name2$Duplicated <- duplicated(df_Name2$Combined)
            df_Name3 <- filter(df_Name2, Duplicated == TRUE)
            df_Name4 <- filter(df_Name2, Combined == df_Name3[1, 5])
            
            if (df_Name4[1,10] > df_Name4[2,10]) {
                
                #Peaks are compared based on UV count first.
                #The peak with the lowest UV count is removed.
                
                Bad_Peak <- df_Name4[2,2]
                
            }   else if (df_Name4[1,10] < df_Name4[2,10]) {
                Bad_Peak <- df_Name4[1,2]
                
            }   else if (df_Name4[1,11] < df_Name4[2,11]) {
                
                #If the UV counts are equal the subtracted UV mean is compared.
                #The peak with the highest UV mean is removed.
                
                Bad_Peak <- df_Name4[2,2]
            }   else {
                Bad_Peak <- df_Name4[1,2]
            }
            df_Name[Bad_Peak, 5:12] <- NA
            df_Name <- select(df_Name, -Combined)
            write.csv(df_Name, 
                      paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                             Coculture_Name, ".CSV"), row.names = FALSE)
            Logic_Row_No <- Logic_Row_No + 1
        }
    }
}

#############################################
#3.NON_UV Peak Matcher: Further assigns peaks based on weaker criteria.
#############################################

Non_UV_Matcher <- function(Interaction_Matrix) {
    
    Matrix_TotalRows <- nrow(Interaction_Matrix)
    Matrix_Row_No <- 1
    
    while (Matrix_Row_No <= Matrix_TotalRows) {
        
        #Reads in the first coculture output file to be amended.
        #Reads in the the corresponding CON files from raw NovaC.
        
        CON1_Name <- as.character(Interaction_Matrix[Matrix_Row_No,1])
        CON2_Name <- as.character(Interaction_Matrix[Matrix_Row_No,2])
        Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
        df_Name <- 
            read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                            Coculture_Name, ".csv"))
        CON1 <- as.data.frame(Read_Excel(CON1_Name))
        CON2 <- as.data.frame(Read_Excel(CON2_Name))
        Coculture <- as.data.frame(Read_Excel(Coculture_Name))
        
        Matrix_Row_No <- Matrix_Row_No + 1
        n <- nrow(df_Name)
        CC_peak <- 1
        
        #'CC_peak' corresponds to the peak no. to be matched in the coculture
        #This is a second set of peak matching that improves on the first round.
        
        while (CC_peak < n+1) {
            CON1_peak <- which(abs(CON1$RetTime-Coculture$RetTime[CC_peak]) ==
                            min(abs(CON1$RetTime-Coculture$RetTime[CC_peak])))
            CON2_peak <- which(abs(CON2$RetTime-Coculture$RetTime[CC_peak]) ==
                            min(abs(CON2$RetTime-Coculture$RetTime[CC_peak])))
            FinalCount1 <- UVcheck(CON1, Coculture, CC_peak, CON1_peak)
            FinalCount2 <- UVCheck(CON2, Coculture, CC_peak, CON2_peak)
            if (!is.na(df_Name$Matched_CON[CC_peak]) | 
                (any(df_Name[,7] == CON1$RetTime[CON1_peak], na.rm = TRUE)) |
                (any(df_Name[,7] == CON2$RetTime[CON2_peak], na.rm = TRUE))) {
                
                #Checks for peak matching already, and skips to the next peak.

            }   else if (Coculture$RetTime[CC_peak] < (CON1$RetTime[CON1_peak] + 0.05) && 
                         Coculture$RetTime[CC_peak] > (CON1$RetTime[CON1_peak] -0.05) &&
                         FinalCount1 < 2)    {
                
                #Checks if the closest match in CON1 satisfies this test.
                
                CON_peak <- CON1_peak
                ratio = (((Coculture[CC_peak,3] - CON1[CON_peak,3])/CON1[CON_peak,3])*100)
                Effect <- Simple_Effect_Categoriser(ratio)
                
                #Assignments of the matched peak
                
                df_Name$Matched_CON[CC_peak] <- CON1_Name
                df_Name$PeakNo_CON[CC_peak] <- CON1$Peak[CON_peak]
                df_Name$RetTime_CON[CC_peak] <- CON1$RetTime[CON_peak]
                df_Name$PeakArea_CON[CC_peak] <- CON1$Area[CON_peak]
                df_Name$UV_Count[CC_peak] <- FinalCount1
                df_Name$PeakRatio[CC_peak] <- ratio
                df_Name$Metabolite_Effect[CC_peak] <- Effect
                
            }   else if (Coculture$RetTime[CC_peak] < (CON2$RetTime[CON2_peak] + 0.05) && 
                         Coculture$RetTime[CC_peak] > (CON2$RetTime[CON2_peak] -0.05) &&
                         FinalCount2 < 2)    {
                CON_peak <- CON2_peak
                ratio = (((Coculture[CC_peak,3] - CON2[CON_peak,3])/CON2[CON_peak,3])*100)
                Effect <- Simple_Effect_Categoriser(ratio)
                df_Name$Matched_CON[CC_peak] <- CON2_Name
                df_Name$PeakNo_CON[CC_peak] <- CON2$Peak[CON_peak]
                df_Name$RetTime_CON[CC_peak] <- CON2$RetTime[CON_peak]
                df_Name$PeakArea_CON[CC_peak] <- CON2$Area[CON_peak]
                df_Name$UV_Count[CC_peak] <- FinalCount2
                df_Name$PeakRatio[CC_peak] <- ratio
                df_Name$Metabolite_Effect[CC_peak] <- Effect
            }    
            CC_peak <- CC_peak + 1
        }
        write.csv(df_Name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         Coculture_Name, ".CSV"), row.names = FALSE)
    }
}

#############################################
#4.Inhibition Checker: Verifies if one culture is inhibited
#############################################

Inhibition_Checker <- function(df_Name, Inhibition_df, Coculture_Name) {
    
    #From an opened output file (might be possible just downstream of missing control peaks)
    match <- as.data.frame(table(df_Name$Matched_CON))
    names(match)[2] <- 'match'
    unmatch <- as.data.frame(table(df_Name$Sample_Ref))
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
                         c("Coculture_Name", "Inhibition", "Dominating_Culture"))
        combined$Coculture_Name <- Coculture_Name
        combined$Var1 <- as.character(combined$Var1)
        temp$Coculture_Name <- combined[1,6]
        temp$Inhibition <- combined[1,5]
        temp$Dominating_Culture <- combined[1,1]
        Inhibition_df <- rbind(Inhibition_df, temp)
    }
    return(Inhibition_df)
}

#############################################
#5.Missing Control Peaks: Adds in the unassigned peaks from the control(s)
#############################################

Missing_Control_Peaks <- function(Interaction_Matrix) {
    
    #Makes a new table that is used in Inhibition_Checker function:
    Inhibition_df <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                              c("Coculture_Name", "Inhibition", "Inhibited_Culture"))
    
    Matrix_TotalRows <- nrow(Interaction_Matrix)
    Matrix_Row_No <- 1
    
    while (Matrix_Row_No <= Matrix_TotalRows) {
        
        #Reads in the first coculture output file to be amended.
        #Reads in the the corresponding CON files from raw NovaC.
        
        CON1_Name <- as.character(Interaction_Matrix[Matrix_Row_No,1])
        CON2_Name <- as.character(Interaction_Matrix[Matrix_Row_No,2])
        Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
        df_Name <- 
            read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                            Coculture_Name, ".csv"))
        
        df_Name$Sample_Ref <- as.character(df_Name$Sample_Ref)
        df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), 
                         sep = "-", remove = FALSE)
        CON1 <- as.data.frame(Read_Excel(CON1_Name))
        CON2 <- as.data.frame(Read_Excel(CON2_Name))
        
        Matrix_Row_No <- Matrix_Row_No + 1
        n <- nrow(CON1)
        CC_peak <- 1
        
        #Sequentially checks CON1 for peak 'CC_peak' in coculture output file
        
        while (CC_peak <= n) {
            
            if (any(df_Name[,5] == paste0(CON1_Name, "-", CC_peak), na.rm = TRUE)) {
                CC_peak <- CC_peak + 1    
            }   else if (any(df_Name[,5] != paste0(CON1_Name, "-", CC_peak), 
                             na.rm = TRUE)) {
                df_Name <- rbind(df_Name, 
                                 c(Sample_Ref = CON1_Name, PeakNo_CC = NA, 
                                   RetTime_CC = NA, PeakArea_CC = NA, 
                                   Combined = NA, MAtched_CON = NA, 
                                   PeakNo_CON = CON1$Peak[CC_peak], 
                                   RetTime_CON = CON1$RetTime[CC_peak], 
                                   PeakArea_CON = CON1$Area[CC_peak], 
                                   UV_Count = NA, Subtracted_UV_Mean = NA, 
                                   PeakRatio = -100, Metabolite_Effect = 1))
            }
            CC_peak <- CC_peak + 1 
        }
        
        n <- nrow(CON2)
        CC_peak <- 1
        
        while (CC_peak <= n) {
            
            if (any(df_Name[,5] == paste0(CON2_Name, "-", CC_peak), na.rm = TRUE)) {
                CC_peak <- CC_peak + 1    
            }   else if (any(df_Name[,5] != paste0(CON2_Name, "-", CC_peak), 
                             na.rm = TRUE)) {
                df_Name <- rbind(df_Name, 
                                 c(Sample_Ref = CON2_Name, PeakNo_CC = NA, 
                                   RetTime_CC = NA, PeakArea_CC = NA,
                                   Combined = NA, MAtched_CON = NA, 
                                   PeakNo_CON = CON2$Peak[CC_peak], 
                                   RetTime_CON = CON2$RetTime[CC_peak], 
                                   PeakArea_CON = CON2$Area[CC_peak], 
                                   UV_Count = NA, Subtracted_UV_Mean = NA, 
                                   PeakRatio = -100, Metabolite_Effect = 1))
                CC_peak <- CC_peak + 1
            }   else {
                CC_peak <- CC_peak + 1
            }
        }
        df_Name <- select(df_Name, -Combined)
        
        #Space to include function to identify inhibition
        Inhibition_df <- Inhibition_Checker(df_Name, Inhibition_df, 
                                            Coculture_Name)
        write.csv(df_Name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         Coculture_Name, ".CSV"), row.names = FALSE)
    }
    Inhibition_df <- transform(Inhibition_df, Inhibition = as.logical(Inhibition))
    write.csv(Inhibition_df, 
              paste0("Testing Broad-Scale Interactions/OutputFiles/Inhibition_df.CSV"), 
              row.names = FALSE)
}

#############################################
#MIMI_Part_2: Carries out the main dereplication processes
#############################################

#Carries out the last two functions:
#1.Double_Peak_Remover: Multiple peaks in a coculture matched to the same
#unique peak of a control
#2.Missing_Control_Peaks: Unique peaks from control(s) not matched to a peak
#in the coculture, are added into a single, unified df

MIMI2 <- function() {  
    
print("Initiating Double_Peak_Remover")
Double_Peak_Remover(Interaction_Matrix)
print("Initiating Non_UV_Matcher")
Non_UV_Matcher(Interaction_Matrix)
print("Initiating Missing_Control_Peaks")
Missing_Control_Peaks(Interaction_Matrix)
print("MIMI_2 completed.")

}
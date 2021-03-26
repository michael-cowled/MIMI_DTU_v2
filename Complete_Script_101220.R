#############################################
#############################################
#Microbial Interaction Metabolite Integrator#
#############################################
#############################################

#A Matrix defined as Interaction_Matrix is required to be made with the THREE
#samples to be compaared.

#R packages required to be loaded in:

library(dplyr)
library(tidyr)
library(readxl)

#Functions to be pre-loaded prior to use of the main function, MIMI()
#1.Read_Excel
#2.UVcheck
#3.UVcheck2
#4.Read_UV
#5.UVSubtract1
#6.UVSubtract2
#7.PeakMatcher1
#8.PeakMatcher2
#9.CON_Consolidator
#10.Effect_Categoriser
#11.Simple_Effect_Categoriser
#12.NON_UV_Peak_Matcher
#13.Double_Peak_Remover
#14.Inhibition_Checker
#15.Missing_Control_Peaks

#Reading in the functions:

#############################################
#1.Read_Excel: Creates 3 data frames for the 3 samples to be compared.
#############################################

Read_Excel <- function(Excel_Name) {
    
    #Input is Excel_Name which is read from the Interaction_Matrix object
    
    #Read in a df based on NovaC excel format
    df_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", 
                                 Excel_Name, ".xlsm"), skip = 3)
    
    #A new df is set up to capture the top 5 UV maxima for a peak
    UV_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
                      c("UV1", "UV2", "UV3", "UV4", "UV5"))
    
    #The following extracts and sorts the top 5 UV maxima from each peak
    n <- nrow(df_Name)
    i <- 1
    while (i <= n) {
        if (!is.na(df_Name$'UV Peaks')) {
            UV_separate <-df_Name$`UV Peaks`[[i]]
            UV_separate <- gsub("\\(","", UV_separate)
            UV_separate <- gsub("\\)","", UV_separate)
            UV_separate <- gsub("< 190","", UV_separate)
            UV_separate <- gsub("s","", UV_separate)
            UV_separate <- strsplit(UV_separate, "\r\n")
            UV_separate <- as.data.frame(UV_separate, col.names = "UV1")
            UV_separate <- separate(UV_separate, 'UV1', c("UV1", "UV1_percent"), 
                                    sep = " ")
            UV_separate <- transform(UV_separate, UV1_percent = as.numeric(UV1_percent))
            UV_separate <- UV_separate[order(UV_separate$UV1_percent, decreasing = TRUE),]
            UV_df2 <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), 
                              c("UV1", "UV2", "UV3", "UV4", "UV5"))
            UV_df2[1,1] <- UV_separate[1,1]
            UV_df2[1,2] <- UV_separate[2,1]
            UV_df2[1,3] <- UV_separate[3,1]
            UV_df2[1,4] <- UV_separate[4,1]
            UV_df2[1,5] <- UV_separate[5,1]
            UV_df <- rbind(UV_df, UV_df2)
            i <- i +1    
        }   else {
            i <- i + 1
            UV_separate <- c(NA, NA, NA, NA, NA)
            UV_df <- rbind(UV_df, UV_separate[1:5,1])
        }
    }
    #UV_df is added onto the loaded df_name object, correlating peaks with UV
    names(UV_df) <- c("UV1", "UV2", "UV3", "UV4", "UV5")   
    UV_df <- transform(UV_df, UV1 = as.numeric(UV1))
    UV_df <- transform(UV_df, UV2 = as.numeric(UV2))
    UV_df <- transform(UV_df, UV3 = as.numeric(UV3))
    UV_df <- transform(UV_df, UV4 = as.numeric(UV4))
    UV_df <- transform(UV_df, UV5 = as.numeric(UV5))
    df_Name <- cbind(df_Name, UV_df)
    df_Name <- select(df_Name, 1:4, 20:24)
    return(df_Name)
}

#############################################
#2/3.UV_Check: Verifies the matching UV maxima for peaks in CON1/2 & Coculture
#############################################

#UVcheck1: Comparing CON1 to coculture

UVcheck1 <- function(CON1, Coculture, i, z) {
    
    #Inputs are CON1 and coculture as read from the Interaction_Matrix
    #'i' refers to the peak no. being compared in Coculture
    #'z' refers to the peak no. being compared in CON1
    
    j <- 5
    k <- 5
    
    #Comparing the 5 UVs from left to right matching to peaks 'i' and 'z'
    #'k' refers to the UV no. being compared in Coculture
    #'j' refers to the UV no. being compared in CON1
    
    uvcount <- 0
    
    #uvcount will be the metric to count the number of matching UVs
    #Checks if UV maxima are within 2 nm of each other
    
    while (k < 10) {
        j<-5
        while (j < 10 && k < 10) {                      #UVs are in columns 5-9
            if (j==9)  {
                j <- j+1
                k <- k+1
            }   else if (is.na(Coculture[i,k])) {
                k <- k+1
            }   else if (is.na(CON1[z,j])) {
                j <- j+1
            }   else if (Coculture[i,k] <= CON1[z,j] + 2 && 
                         Coculture[i,k] >= CON1[z,j] - 2) {
                uvcount <- uvcount + 1
                j <- j+1
            }   else {
                j <- j+1
            }
        }
    }
    return(uvcount)
}

#UVcheck2: Comparing CON2 to coculture

UVCheck2 <- function(CON2, Coculture, i, z) {
    
    #Inputs are CON2 and coculture as read from the Interaction_Matrix
    #'i' refers to the peak no. being compared in Coculture
    #'z' refers to the peak no. being compared in CON2
    
    uvcount <- 0
    j <- 5
    k <- 5 
    while (k < 10) {
        j<-5
        while (j < 10 && k < 10) {
            if (j==9)  {
                j <- j+1
                k <- k+1
            }   else if (is.na(Coculture[i,k])) {
                k <- k+1
            }   else if (is.na(CON2[z,j])) {
                j <- j+1
            }   else if (Coculture[i,k] <= CON2[z,j] + 2 && 
                         Coculture[i,k] >= CON2[z,j] - 2) {
                uvcount <- uvcount + 1
                j <- j+1
            }   else {
                j <- j+1
            }
        }
    }
    return(uvcount)
}


#############################################
#4.Read_UV Function: Reads in the raw UV data for every wavelength
#############################################

Read_UV <- function(Excel_Name) {
    
    #Input is Excel_Name which is read from the Interaction_Matrix object
    #Generates a df with col1 = Wavelength, col2 = Abs for Peak 1, etc.
    
    UV_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/",
                                 Excel_Name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(UV_Name)
}

#############################################
#5/6.UV Subtract: Subtracts whole UV spectra and calculates the mean for peaks 
#in CON1/2 & Coculture
#############################################

#UVsubtract1: Comparing CON1 to coculture

UVSubtract1 <- function(CON1_UV, Coculture_UV, i, z) {
    
    #Inputs are CON1_UV and coculture_UV as created from the Interaction_Matrix
    #'i' refers to the peak no. being compared in Coculture
    #'z' refers to the peak no. being compared in CON2
    
    k <- i + 1
    j <- z + 1
    
    #Subtracts absorbances at each wavelength and calculates the column mean
    #'k' refers to the column of abs data being compared in Coculture
    #'j' refers to the column of abs data being compared in CON1
    #Note: Col2 = Peak#1, Col3 = Peak#2, etc. hence '+1'
    
    SubtractedUV <- Coculture_UV[,k]-CON1_UV[,j]
    UV_Mean <- mean(SubtractedUV)
    return(UV_Mean)
}

#UVsubtract2: Comparing CON2 to coculture

UVSubtract2 <- function(CON2_UV, Coculture_UV, i, z) {
    k <- i + 1
    j <- z + 1
    SubtractedUV <- Coculture_UV[,k]-CON2_UV[,j]
    UV_Mean <- mean(SubtractedUV)
    return(UV_Mean)
}

#############################################
#7.Peak Matcher: Finds and compares the nearest matching peak in CON1
#############################################

Peak_Matcher1 <- function(CON1, Coculture, CON1_UV, Coculture_UV) {
    
    Coculture_df <- data.frame(PeakNo_CC = NA, RetTime_CC = NA, 
                               PeakArea_CC = NA, PeakNo_CON1 = NA, 
                               RetTime_CON1 = NA, PeakArea_CON1 = NA,
                               UV_Count_CON1 = NA, 
                               Subtracted_UV_Mean_CON1 = NA, 
                               PeakRatio_CON1 = NA)
    
    n <- nrow(Coculture)
    i = 1
    
    #'i' corresponds to the peak no. to be compared in the coculture
    
    while (i < n+1) {
        z <- which(abs(CON1$RetTime-Coculture$RetTime[i]) ==
                       min(abs(CON1$RetTime-Coculture$RetTime[i])))
        
        #'z' finds the peak in CON1 with the closest ret-time to peak 'i'
        
        ratio = (((Coculture[i,3] - CON1[z,3])/CON1[z,3])*100) 
        
        #Computes the ratio of peak areas as a %
        
        FinalCount <- UVcheck1(CON1, Coculture, i, z)
        UV_Mean <- UVSubtract1(CON1_UV, Coculture_UV, i, z)
        
        #Performs both the UVcheck and UVsubtract functions
        
        if (Coculture$RetTime[i] < (CON1$RetTime[z] + 0.15) && 
            Coculture$RetTime[i] > (CON1$RetTime[z] -0.15) &&
            (FinalCount > 2 || abs(UV_Mean) < 1.5)) {
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UVcount is at least 3 OR UV_Mean < 1.5
            #If satisfied, the peak is declared a match and added to df
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[i], 
                                    RetTime_CC = Coculture$RetTime[i], 
                                    PeakArea_CC = Coculture$Area[i], 
                                    PeakNo_CON1 = CON1$Peak[z], 
                                    RetTime_CON1 = CON1$RetTime[z],
                                    PeakArea_CON1 = CON1$Area[z], 
                                    UV_Count_CON1 = FinalCount, 
                                    Subtracted_UV_Mean_CON1 = UV_Mean, 
                                    PeakRatio_CON1 = ratio))
            i <- i+1
        }    else if (Coculture$RetTime[i] < (CON1$RetTime[z] + 0.15) && 
                      Coculture$RetTime[i] > (CON1$RetTime[z] -0.15) &&
                      FinalCount > 1 && abs(UV_Mean) < 2) {  
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UVcount is at least 2 AND UV_Mean < 2
            #Both UVcount and UV_Mean are less strict requirement
            #But requiring satisfaction of both adds stringency
            #If satisfied, the peak is declared a match and added to df
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[i], 
                                    RetTime_CC = Coculture$RetTime[i], 
                                    PeakArea_CC = Coculture$Area[i], 
                                    PeakNo_CON1 = CON1$Peak[z], 
                                    RetTime_CON1 = CON1$RetTime[z],
                                    PeakArea_CON1 = CON1$Area[z], 
                                    UV_Count_CON1 = FinalCount, 
                                    Subtracted_UV_Mean_CON1 = UV_Mean, 
                                    PeakRatio_CON1 = ratio))
            i <- i+1
            
        }    else if (Coculture$RetTime[i] < (CON1$RetTime[z] + 0.05) && 
                      Coculture$RetTime[i] > (CON1$RetTime[z] -0.05) &&
                      FinalCount > 1) {  
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[i], 
                                    RetTime_CC = Coculture$RetTime[i], 
                                    PeakArea_CC = Coculture$Area[i], 
                                    PeakNo_CON1 = CON1$Peak[z], 
                                    RetTime_CON1 = CON1$RetTime[z],
                                    PeakArea_CON1 = CON1$Area[z], 
                                    UV_Count_CON1 = FinalCount, 
                                    Subtracted_UV_Mean_CON1 = UV_Mean, 
                                    PeakRatio_CON1 = ratio))
            i <- i+1
            
        }   else  { 
            
            #The peak is declared NOT to have a match in CON1, and NA added
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[i], 
                                    RetTime_CC = Coculture$RetTime[i], 
                                    PeakArea_CC = Coculture$Area[i], 
                                    PeakNo_CON1 = NA, RetTime_CON1 = NA, 
                                    PeakArea_CON1 = NA, UV_Count_CON1 = NA, 
                                    Subtracted_UV_Mean_CON1 = NA, 
                                    PeakRatio_CON1 = NA))
            i <- i+1
        }   
    }
    return(Coculture_df)
}

#############################################
#8.Peak Matcher2: Finds and compares the nearest matching peak in CON2
#############################################

Peak_Matcher2 <- function(CON2, Coculture, CON2_UV, Coculture_UV) {
    
    Coculture_df2 <- data.frame(PeakNo_CC = NA, RetTime_CC = NA,
                                PeakArea_CC = NA, PeakNo_CON2 = NA, 
                                RetTime_CON2 = NA, PeakArea_CON2 = NA,
                                UV_Count_CON2 = NA, 
                                Subtracted_UV_Mean_CON2 = NA, 
                                PeakRatio_CON2 = NA)
    n <- nrow(Coculture)
    i = 1
    
    while (i < n+1) {
        z <- which(abs(CON2$RetTime-Coculture$RetTime[i]) ==
                       min(abs(CON2$RetTime-Coculture$RetTime[i])))
        ratio = (((Coculture[i,3] - CON2[z,3])/CON2[z,3])*100)
        FinalCount <- UVCheck2(CON2, Coculture, i, z)
        UV_Mean <- UVSubtract2(CON2_UV, Coculture_UV, i, z)
        
        if (Coculture$RetTime[i] < (CON2$RetTime[z] + 0.15) && 
            Coculture$RetTime[i] > (CON2$RetTime[z] -0.15)
            && (FinalCount > 2 || abs(UV_Mean) < 1.5)) {  
            Coculture_df2 <- rbind(Coculture_df2, 
                                   c(PeakNo_CC = Coculture$Peak[i], 
                                     RetTime_CC = Coculture$RetTime[i], 
                                     PeakArea_CC = Coculture$Area[i], 
                                     PeakNo_CON2 = CON2$Peak[z], 
                                     RetTime_CON2 = CON2$RetTime[z], 
                                     PeakArea_CON2 = CON2$Area[z], 
                                     UV_Count_CON2 = FinalCount, 
                                     Subtracted_UV_Mean_CON2 = UV_Mean, 
                                     PeakRatio_CON2 = ratio))
            i <- i+1
        }   else if (Coculture$RetTime[i] < (CON2$RetTime[z] + 0.15) && 
                     Coculture$RetTime[i] > (CON2$RetTime[z] -0.15) &&
                     FinalCount > 1 && abs(UV_Mean) < 2) {  
            Coculture_df2 <- rbind(Coculture_df2, 
                                   c(PeakNo_CC = Coculture$Peak[i], 
                                     RetTime_CC = Coculture$RetTime[i], 
                                     PeakArea_CC = Coculture$Area[i], 
                                     PeakNo_CON2 = CON2$Peak[z], 
                                     RetTime_CON2 = CON2$RetTime[z], 
                                     PeakArea_CON2 = CON2$Area[z], 
                                     UV_Count_CON2 = FinalCount, 
                                     Subtracted_UV_Mean_CON2 = UV_Mean, 
                                     PeakRatio_CON2 = ratio))
            i <- i+1
        }   else if (Coculture$RetTime[i] < (CON2$RetTime[z] + 0.05) && 
                     Coculture$RetTime[i] > (CON2$RetTime[z] -0.05) &&
                     FinalCount > 1) {  
            Coculture_df2 <- rbind(Coculture_df2, 
                                   c(PeakNo_CC = Coculture$Peak[i], 
                                     RetTime_CC = Coculture$RetTime[i], 
                                     PeakArea_CC = Coculture$Area[i], 
                                     PeakNo_CON2 = CON2$Peak[z], 
                                     RetTime_CON2 = CON2$RetTime[z], 
                                     PeakArea_CON2 = CON2$Area[z], 
                                     UV_Count_CON2 = FinalCount, 
                                     Subtracted_UV_Mean_CON2 = UV_Mean, 
                                     PeakRatio_CON2 = ratio))
            i <- i+1
        }   else  { 
            Coculture_df2 <- rbind(Coculture_df2, 
                                   c(PeakNo_CC = Coculture$Peak[i], 
                                     RetTime_CC = Coculture$RetTime[i], 
                                     PeakArea_CC = Coculture$Area[i], 
                                     PeakNo_CON2 = NA, RetTime_CON2 = NA, 
                                     PeakArea_CON2 = NA, UV_Count_CON2 = NA,
                                     Subtracted_UV_Mean_CON2 = NA, 
                                     PeakRatio_CON2 = NA))
            i <- i+1
        }   
    }
    return(Coculture_df2)
}

#############################################
#9.CON Consolidator: Removes double peak matching to a coculture peak
#############################################

#The following code removes double assignments to a coculture peak
#from two different controls.

#An example is that perhaps peak 1 from the coculture matches
#to peak 1 in CON1 and peak 2 in CON2.

#This function will determine which control's matched peak matches
#closest and remove the assignment for the weakest match.

CON_Consolidator <- function(Coculture_df, CON1_Name, CON2_Name) {
    
    RowNo <- 1
    Rows <- nrow(Coculture_df)
    
    #sets up a new tidier version of the df
    
    MatchedPeak_df <- setNames(data.frame(matrix(ncol = 7, 
                                                 nrow = nrow(Coculture_df))), 
                               c("Matched_CON", "PeakNo_CON", "RetTime_CON", 
                                 "PeakArea_CON", "UV_Count", 
                                 "Subtracted_UV_Mean", "PeakRatio"))
    
    while (RowNo <= Rows+1) {
        if (is.na(Coculture_df[RowNo, 4]) && 
            is.na(Coculture_df[RowNo, 10])) {
            RowNo <- RowNo + 1
        }   else if (!is.na(Coculture_df[RowNo, 4]) && 
                     !is.na(Coculture_df[RowNo, 10] && 
                            abs(Coculture_df[RowNo, 8]) < 
                            abs(Coculture_df[RowNo, 14]))) { 
            
            #When there are two peaks matched (!is.na for both),
            #the UV_means are compared, with the higher removed.
            
            MatchedPeak_df[RowNo, 1] <- CON1_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
            RowNo <- RowNo + 1
        }   else if (!is.na(Coculture_df[RowNo, 4]) && 
                     !is.na(Coculture_df[RowNo, 10] && 
                            abs(Coculture_df[RowNo, 8]) > 
                            abs(Coculture_df[RowNo, 14]))) { 
            MatchedPeak_df[RowNo, 1] <- CON2_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
            RowNo <- RowNo + 1
        }   else if (is.na(Coculture_df[RowNo, 10]))   {
            
            #If no double-peak mactching but signifies a matched CON1 peak
            #then CON1 peak set as the matched peak
            
            MatchedPeak_df[RowNo, 1] <- CON1_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
            RowNo <- RowNo + 1
        }   else {
            
            #If no double-peak mactching or matched CON1 peak
            #then CON2 peak set as the matched peak
            
            MatchedPeak_df[RowNo, 1] <- CON2_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
            RowNo <- RowNo + 1
        }
    }
    return(MatchedPeak_df)
}

#############################################
#10.Metabolite Effect Characteriser: Converting PeakRatio to a Factor
#############################################

#The following piece of code adds a column that categorises the peak 
#areas into suppressions, and enhancements

Effect_Categoriser <- function(Refined_Coculture_df, Coculture_Name) {
    
    RowNo <- 1
    Rows <- nrow(Refined_Coculture_df)
    
    #A df is created to list the effects corresponding to matched peaks
    
    Metabolite_effect_df <- setNames(data.frame(matrix(ncol=1, nrow=Rows)), 
                                     c("Metabolite_Effect"))
    
    while (RowNo <= Rows) {
        if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
            is.na(Refined_Coculture_df[RowNo, 11])) {
            Metabolite_effect_df[RowNo, 1] <- 6 #Induction
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] > -100
                     && Refined_Coculture_df[RowNo, 11] <= -20) {
            Metabolite_effect_df[RowNo, 1] <- 2 #Suppression
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] > -20
                     && Refined_Coculture_df[RowNo, 11] < 20) {
            Metabolite_effect_df[RowNo, 1] <- 3 #Little to No Change
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] >= 20
                     && Refined_Coculture_df[RowNo, 11] < 100) {
            Metabolite_effect_df[RowNo, 1] <- 4 #Enhancement
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] >= 100) {
            Metabolite_effect_df[RowNo, 1] <- 5 #Major Enhancement
            RowNo <- RowNo + 1
        }
    }
    Refined_Coculture_df <- 
        cbind(Refined_Coculture_df, Metabolite_effect_df)
    return(Refined_Coculture_df)
}

#############################################
#11.Simple Metabolite Effect Characteriser:
#############################################

#This is a simple version of the effect categoriser with a single ratio input

Simple_Effect_Categoriser <- function(ratio) {
    if (ratio > -100 && ratio <= 20) {
        Effect <- 2
        return(Effect)
    }   else if (ratio > -20 && ratio < 20) {
        Effect <- 3
        return(Effect)
    }   else if (ratio >= 20 && ratio < 100) {
        Effect <- 4
        return(Effect)
    }   else if (ratio >= 100) {
        Effect <- 5
        return(Effect)
    }
}

#############################################
#12.NON_UV Peak Matcher: Further assigns peaks based on weaker criteria.
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
        i <- 1
        
        #'i' corresponds to the peak no. to be matched in the coculture
        #This is a second set of peak matching that improves on the first round.
        
        while (i < n+1) {
            z1 <- which(abs(CON1$RetTime-Coculture$RetTime[i]) ==
                            min(abs(CON1$RetTime-Coculture$RetTime[i])))
            z2 <- which(abs(CON2$RetTime-Coculture$RetTime[i]) ==
                            min(abs(CON2$RetTime-Coculture$RetTime[i])))
            FinalCount1 <- UVcheck1(CON1, Coculture, i, z = z1)
            FinalCount2 <- UVCheck2(CON2, Coculture, i, z = z2)
            if (!is.na(df_Name$Matched_CON[i]) | 
                (any(df_Name[,7] == CON1$RetTime[z1], na.rm = TRUE)) |
                (any(df_Name[,7] == CON2$RetTime[z2], na.rm = TRUE))) {
                
                #Checks for peak matching already, and skips to the next peak.
                
                i <- i+1
            }   else if (Coculture$RetTime[i] < (CON1$RetTime[z1] + 0.05) && 
                         Coculture$RetTime[i] > (CON1$RetTime[z1] -0.05) &&
                         FinalCount1 < 2)    {
                
                #Checks if the closest match in CON1 satisfies this test.
                
                z <- z1
                ratio = (((Coculture[i,3] - CON1[z,3])/CON1[z,3])*100)
                Effect <- Simple_Effect_Categoriser(ratio)
                
                #Assignments of the matched peak
                
                df_Name$Matched_CON[i] <- CON1_Name
                df_Name$PeakNo_CON[i] <- CON1$Peak[z]
                df_Name$RetTime_CON[i] <- CON1$RetTime[z]
                df_Name$PeakArea_CON[i] <- CON1$Area[z]
                df_Name$UV_Count[i] <- FinalCount1
                df_Name$PeakRatio[i] <- ratio
                df_Name$Metabolite_Effect[i] <- Effect
                
                i <- i+1
                
            }   else if (Coculture$RetTime[i] < (CON2$RetTime[z2] + 0.05) && 
                         Coculture$RetTime[i] > (CON2$RetTime[z2] -0.05) &&
                         FinalCount2 < 2)    {
                z <- z2
                ratio = (((Coculture[i,3] - CON2[z,3])/CON2[z,3])*100)
                Effect <- Simple_Effect_Categoriser(ratio)
                df_Name$Matched_CON[i] <- CON2_Name
                df_Name$PeakNo_CON[i] <- CON2$Peak[z]
                df_Name$RetTime_CON[i] <- CON2$RetTime[z]
                df_Name$PeakArea_CON[i] <- CON2$Area[z]
                df_Name$UV_Count[i] <- FinalCount2
                df_Name$PeakRatio[i] <- ratio
                df_Name$Metabolite_Effect[i] <- Effect
                i <- i+1
            }    else {  
                i <- i+1
            }
        }
        write.csv(df_Name, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         Coculture_Name, ".CSV"), row.names = FALSE)
    }
}

#############################################
#13.Double Peak Remover: Removes doubley-assigned peaks from a control.
#############################################

#An example is that perhaps the peaks 1 and 2 from the coculture match to
#twice to peak 1 in the control.

#This function will determine which peak from the control matches closest
#and remove the assignment for the weakest match.

Double_Peak_Remover <- function(Interaction_Matrix) {
    
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
            df_Name4 <- df_Name
            df_Name <- filter(df_Name, Combined != "NA-NA")
            df_Name$Duplicated <- duplicated(df_Name$Combined)
            df_Name2 <- filter(df_Name, Duplicated == TRUE)
            df_Name3 <- filter(df_Name, Combined == df_Name2[1, 5])
            
            if (df_Name3[1,10] > df_Name3[2,10]) {
                
                #Peaks are compared based on UV count first.
                #The peak with the lowest UV count is removed.
                
                Bad_Peak <- df_Name3[2,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }   else if (df_Name3[1,10] < df_Name3[2,10]) {
                Bad_Peak <- df_Name3[1,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }   else if (df_Name3[1,11] < df_Name3[2,11]) {
                
                #If the UV counts are equal the subtracted UV mean is compared.
                #The peak with the highest UV mean is removed.
                
                Bad_Peak <- df_Name3[2,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }   else {
                Bad_Peak <- df_Name3[1,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }
            df_Name4 <- select(df_Name4, -Combined)
            write.csv(df_Name4, 
                      paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                             Coculture_Name, ".CSV"), row.names = FALSE)
            Logic_Row_No <- Logic_Row_No + 1
        }
    }
}

#############################################
#14.Inhibition Checker: Verifies if one culture is inhibited
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
#15.Missing Control Peaks: Adds in the unassigned peaks from the control(s)
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
        i <- 1
        
        #Sequentially checks CON1 for peak 'i' in coculture output file
        
        while (i <= n) {
            
            if (any(df_Name[,5] == paste0(CON1_Name, "-", i), na.rm = TRUE)) {
                i <- i +1    
            }   else if (any(df_Name[,5] != paste0(CON1_Name, "-", i), 
                             na.rm = TRUE)) {
                df_Name <- rbind(df_Name, 
                                 c(Sample_Ref = CON1_Name, PeakNo_CC = NA, 
                                   RetTime_CC = NA, PeakArea_CC = NA, 
                                   Combined = NA, MAtched_CON = NA, 
                                   PeakNo_CON = CON1$Peak[i], 
                                   RetTime_CON = CON1$RetTime[i], 
                                   PeakArea_CON = CON1$Area[i], 
                                   UV_Count = NA, Subtracted_UV_Mean = NA, 
                                   PeakRatio = -100, Metabolite_Effect = 1))
                i <- i + 1
            }   else {
                i <- i + 1
            }
        }
        
        n <- nrow(CON2)
        i <- 1
        
        while (i <= n) {
            
            if (any(df_Name[,5] == paste0(CON2_Name, "-", i), na.rm = TRUE)) {
                i <- i +1    
            }   else if (any(df_Name[,5] != paste0(CON2_Name, "-", i), 
                             na.rm = TRUE)) {
                df_Name <- rbind(df_Name, 
                                 c(Sample_Ref = CON2_Name, PeakNo_CC = NA, 
                                   RetTime_CC = NA, PeakArea_CC = NA,
                                   Combined = NA, MAtched_CON = NA, 
                                   PeakNo_CON = CON2$Peak[i], 
                                   RetTime_CON = CON2$RetTime[i], 
                                   PeakArea_CON = CON2$Area[i], 
                                   UV_Count = NA, Subtracted_UV_Mean = NA, 
                                   PeakRatio = -100, Metabolite_Effect = 1))
                i <- i + 1
            }   else {
                i <- i + 1
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
    
    write.csv(Inhibition_df, 
        paste0("Testing Broad-Scale Interactions/OutputFiles/Inhibition_df.CSV"), 
    row.names = FALSE)
}

#############################################
#MIMI: The main working-function to compare the peak-matching and refinement
#############################################

MIMI <- function() {           
    
    #Will not work unless Interaction_Matrix object created by the user.
    
    Matrix_TotalRows <- nrow(Interaction_Matrix)
    Matrix_Row_No <- 1
    
    while (Matrix_Row_No <= Matrix_TotalRows) {
        
        CON1_Name <- as.character(Interaction_Matrix[Matrix_Row_No,1])
        CON2_Name <- as.character(Interaction_Matrix[Matrix_Row_No,2])
        Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
        Matrix_Row_No <- Matrix_Row_No +1
        
        #Generates 3 dataframes using the Read_Excel function for the first
        #interction to be investigated from the Interaction_Matrix
        
        CON1 <- as.data.frame(Read_Excel(CON1_Name))
        CON2 <- as.data.frame(Read_Excel(CON2_Name))
        Coculture <- as.data.frame(Read_Excel(Coculture_Name))
        
        #Generates 3 dataframes using the Read_UV function
        
        CON1_UV <- as.data.frame(Read_UV(CON1_Name))
        CON2_UV <- as.data.frame(Read_UV(CON2_Name))
        Coculture_UV <- as.data.frame(Read_UV(Coculture_Name))
        
        #sets up a df to be expanded and exported to csv
        
        Coculture_df <- Peak_Matcher1(CON1, Coculture, CON1_UV, Coculture_UV)
        
        #The following repeats for CON2
        
        Coculture_df2 <- Peak_Matcher2(CON2, Coculture, CON2_UV, Coculture_UV)
        
        #The two dfs for CON1 and CON2 are merged together
        
        Coculture_df <- cbind(Coculture_df, Coculture_df2[, 4:9])
        
        #Removes the first arbitrary row of missing values
        
        Coculture_df <- Coculture_df[2:nrow(Coculture_df), ]
        
        #Performs a function to correct for double peak matching to a unique
        #peak to peaks from more than one CON
        
        MatchedPeak_df <- CON_Consolidator(Coculture_df, CON1_Name, CON2_Name)
        
        #Final processing steps in creating the tidied df (Refined_Coculture_df)
        
        MatchedPeak_df <- MatchedPeak_df[1:nrow(MatchedPeak_df), ]
        df <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(Coculture_df))), 
                       "Sample_Ref")
        df[1:nrow(df), ] <- Coculture_Name
        Refined_Coculture_df <- Coculture_df[1:nrow(Coculture_df), 1:3]
        Refined_Coculture_df <- cbind(df, Refined_Coculture_df)
        Refined_Coculture_df <- cbind(Refined_Coculture_df, MatchedPeak_df)
        
        #Performs a function to characterise effects onto metabolites
        #and adds this into the tidied data set
        
        Refined_Coculture_df <- 
            Effect_Categoriser(Refined_Coculture_df, Coculture_Name)
        
        #Rewrites the tidied dataset to file
        
        write.csv(Refined_Coculture_df, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         Coculture_Name, ".CSV"), row.names = FALSE)
    }
    
    #Carries out the last two functions:
    #1.Double_Peak_Remover: Multiple peaks in a coculture matched to the same
    #unique peak of a control
    #2.Missing_Control_Peaks: Unique peaks from control(s) not matched to a peak
    #in the coculture, are added into a single, unified df
    
    print("Initiating Double_Peak_Remover")
    Double_Peak_Remover(Interaction_Matrix)
    print("Initiating Non_UV_Matcher")
    Non_UV_Matcher(Interaction_Matrix)
    print("Initiating Missing_Control_Peaks")
    Missing_Control_Peaks(Interaction_Matrix)
    print("MIMI completed.")
}

#Load in excel file named Interaction_Matrix and use MIMI() function:
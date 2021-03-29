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
#1.Read_Excel
#2.UVcheck
#3.Read_UV
#4.UVSubtract
#5.PeakMatcher
#6.CON_Consolidator
#7.Effect_Categoriser

#Reading in the functions:

#############################################
#1.Read_Excel: Creates 3 data frames for the 3 samples to be compared.
#############################################

Read_Excel <- function(Excel_Name) {
    
    #Input is Excel_Name which is read from the Interaction_Matrix object
    
    #Read in a df based on NovaC excel format
    raw_df <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", 
                                 Excel_Name, ".xlsm"), skip = 3)
    
    #A new df is set up to capture the top 5 UV maxima for each peak
    UV_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
                      c("UV1", "UV2", "UV3", "UV4", "UV5"))
    
    #The following extracts and sorts the top 5 UV maxima from each peak
    n <- nrow(raw_df)
    CC_peak<- 1
    while (CC_peak<= n) {
        if (!is.na(raw_df$'UV Peaks')) {
            raw_UV <-raw_df$`UV Peaks`[[CC_peak]]
            raw_UV_processed <- gsub("\\(","", raw_UV)
            raw_UV_processed <- gsub("\\)","", raw_UV_processed)
            raw_UV_processed <- gsub("< 190","", raw_UV_processed)
            raw_UV_processed <- gsub("s","", raw_UV_processed)
            raw_UV_split <- strsplit(raw_UV_processed, "\r\n")
            raw_UV_split <- as.data.frame(raw_UV_split, col.names = "UV1")
            raw_UV_split <- separate(raw_UV_split, 'UV1', c("UV1", "UV1_percent"), 
                                    sep = " ")
            raw_UV_split <- transform(raw_UV_split, UV1_percent = as.numeric(UV1_percent))
            raw_UV_ordered <- raw_UV_split[order(raw_UV_split$UV1_percent, decreasing = TRUE),]
            
            #A second UV_df is set up to capture the UV data for a single peak 
            UV_df_peak <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), 
                               c("UV1", "UV2", "UV3", "UV4", "UV5"))
            
            #Transposes UVs into row format; more efficient method possible
            UV_df_peak[1,1] <- raw_UV_ordered[1,1]
            UV_df_peak[1,2] <- raw_UV_ordered[2,1]
            UV_df_peak[1,3] <- raw_UV_ordered[3,1]
            UV_df_peak[1,4] <- raw_UV_ordered[4,1]
            UV_df_peak[1,5] <- raw_UV_ordered[5,1]
            
            #Lists the UV for the associated peak into the originally setup UV_df
            UV_df <- rbind(UV_df, UV_df_peak)
            
        }   else {
            
            #Lists the UVs as nothing instead
            UV_df_peak <- c(NA, NA, NA, NA, NA)
            UV_df <- rbind(UV_df, UV_df_peak)
            
        }   
        CC_peak<- CC_peak+ 1
    }
    
    #UV_df is added onto the loaded raw_df object, correlating peaks with UV
    names(UV_df) <- c("UV1", "UV2", "UV3", "UV4", "UV5")   
    UV_df <- transform(UV_df, UV1 = as.numeric(UV1))
    UV_df <- transform(UV_df, UV2 = as.numeric(UV2))
    UV_df <- transform(UV_df, UV3 = as.numeric(UV3))
    UV_df <- transform(UV_df, UV4 = as.numeric(UV4))
    UV_df <- transform(UV_df, UV5 = as.numeric(UV5))
    raw_df_with_sorted_UV <- cbind(raw_df, UV_df)
    raw_df_with_sorted_UV <- select(raw_df_with_sorted_UV, 1:4, 20:24)
    return(raw_df_with_sorted_UV)
    
}

#############################################
#2.UV_Check: Verifies the matching UV maxima for peaks in CON1/2 & Coculture
#############################################

#UVcheck: Comparing CON1 to coculture
#Inputs are CON1/CON2 and coculture as read from the Interaction_Matrix
#'CC_peak' refers to the peak no. being compared in Coculture
#'CON_peak' refers to the peak no. being compared in control
#UV data is located in columns 5:9, hence 'CON_UV_no' and 'CC_UV_no' set to start at 5.

UVCheck <- function(control, Coculture, CC_peak, CON_peak) {
   
    uvcount <- 0
    CC_UV_no <- 5 
    CON_UV_no <- 5
    
    while (CC_UV_no < 10) {
        CON_UV_no <- 5
        while (CON_UV_no < 10 && CC_UV_no < 10) {
            if (CON_UV_no == 9)  {
                CON_UV_no <- CON_UV_no + 1
                CC_UV_no <- CC_UV_no + 1
            }   else if (is.na(Coculture[CC_peak,CC_UV_no])) {
                CC_UV_no <- CC_UV_no + 1
            }   else if (is.na(control[CON_peak,CON_UV_no])) {
                CON_UV_no <- CON_UV_no + 1
            }   else if (Coculture[CC_peak,CC_UV_no] <= control[CON_peak,CON_UV_no] + 2 && 
                         Coculture[CC_peak,CC_UV_no] >= control[CON_peak,CON_UV_no] - 2) {
                uvcount <- uvcount + 1
                CON_UV_no <- CON_UV_no + 1
            }   else {
                CON_UV_no <- CON_UV_no + 1
            }
        }
    }
    return(uvcount)
}

#############################################
#3.Read_UV Function: Reads in the raw UV data for every wavelength
#############################################

#Input is Excel_Name which is read from the Interaction_Matrix object
#Generates a df with col1 = Wavelength, col2 = Abs for Peak 1, etc.

Read_UV <- function(Excel_Name) {
    UV_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/",
                                 Excel_Name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(UV_Name)
}

#############################################
#4.UV Subtract: Subtracts whole UV spectra and calculates the mean for peaks 
#in CON1/2 & Coculture
#############################################

#UVsubtract: Comparing control to coculture

UVSubtract <- function(CON_UV, Coculture_UV, CC_peak, CON_peak) {
    
    #Inputs are CON1_UV and coculture_UV as created from the Interaction_Matrix
    #'CC_peak' refers to the peak no. being compared in Coculture
    #'CON_peak' refers to the peak no. being compared in CONtrol
    
    CC_UV_no <- CC_peak + 1
    CON_UV_no <- CON_peak + 1
    
    #Subtracts absorbances at each wavelength and calculates the column mean
    #'CC_UV_no' refers to the column of abs data being compared in Coculture
    #'CON_UV_no' refers to the column of abs data being compared in control
    #Note: Col2 = Peak#1, Col3 = Peak#2, etc. hence '+1'
    
    SubtractedUV <- Coculture_UV[,CC_UV_no]-CON_UV[,CON_UV_no]
    UV_Mean <- mean(SubtractedUV)
    return(UV_Mean)
}

#############################################
#5.Peak Matcher: Finds and compares the nearest matching peak in control
#############################################

Peak_Matcher <- function(CON, Coculture, CON_UV, Coculture_UV) {
    
    #Sets up a df to with the desired column names.
    Coculture_df <- data.frame(PeakNo_CC = NA, RetTime_CC = NA, 
                               PeakArea_CC = NA, PeakNo_CON = NA, 
                               RetTime_CON = NA, PeakArea_CON = NA,
                               UV_Count_CON = NA, 
                               Subtracted_UV_Mean_CON = NA, 
                               PeakRatio_CON = NA)
    
    n <- nrow(Coculture)
    CC_peak <- 1
    
    #'CC_peak' corresponds to the peak no. to be compared in the coculture
    
    while (CC_peak< n + 1) {
        CON_peak <- which(abs(CON$RetTime-Coculture$RetTime[CC_peak]) ==
                       min(abs(CON$RetTime-Coculture$RetTime[CC_peak])))
        
        #'CON_peak' finds the peak in CON1 with the closest ret-time to peak 'CC_peak'
        
        ratio = (((Coculture[CC_peak,3] - CON[CON_peak,3])/CON[CON_peak,3])*100) 
        
        #Computes the ratio of peak areas as a %
        
        FinalCount <- UVcheck(CON, Coculture, CC_peak, CON_peak)
        UV_Mean <- UVSubtract(CON_UV, Coculture_UV, CC_peak, CON_peak)
        
        #Performs both the UVcheck and UVsubtract functions
        
        if (Coculture$RetTime[CC_peak] < (CON$RetTime[CON_peak] + 0.15) && 
            Coculture$RetTime[CC_peak] > (CON$RetTime[CON_peak] -0.15) &&
            (FinalCount > 2 || abs(UV_Mean) < 1.5)) {
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UVcount is at least 3 OR UV_Mean < 1.5
            #If satisfied, the peak is declared a match and added to df
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[CC_peak], 
                                    RetTime_CC = Coculture$RetTime[CC_peak], 
                                    PeakArea_CC = Coculture$Area[CC_peak], 
                                    PeakNo_CON = CON$Peak[CON_peak], 
                                    RetTime_CON = CON$RetTime[CON_peak],
                                    PeakArea_CON = CON$Area[CON_peak], 
                                    UV_Count_CON = FinalCount, 
                                    Subtracted_UV_Mean_CON = UV_Mean, 
                                    PeakRatio_CON = ratio))
            
        }    else if (Coculture$RetTime[CC_peak] < (CON$RetTime[CON_peak] + 0.15) && 
                      Coculture$RetTime[CC_peak] > (CON$RetTime[CON_peak] -0.15) &&
                      FinalCount > 1 && abs(UV_Mean) < 2) {  
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UVcount is at least 2 AND UV_Mean < 2
            #Both UVcount and UV_Mean are less strict requirement
            #But requiring satisfaction of both adds stringency
            #If satisfied, the peak is declared a match and added to df
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[CC_peak], 
                                    RetTime_CC = Coculture$RetTime[CC_peak], 
                                    PeakArea_CC = Coculture$Area[CC_peak], 
                                    PeakNo_CON = CON$Peak[CON_peak], 
                                    RetTime_CON = CON$RetTime[CON_peak],
                                    PeakArea_CON = CON$Area[CON_peak], 
                                    UV_Count_CON = FinalCount, 
                                    Subtracted_UV_Mean_CON = UV_Mean, 
                                    PeakRatio_CON = ratio))
            
        }    else if (Coculture$RetTime[CC_peak] < (CON$RetTime[CON_peak] + 0.05) && 
                      Coculture$RetTime[CC_peak] > (CON$RetTime[CON_peak] -0.05) &&
                      FinalCount > 1) {  
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[CC_peak], 
                                    RetTime_CC = Coculture$RetTime[CC_peak], 
                                    PeakArea_CC = Coculture$Area[CC_peak], 
                                    PeakNo_CON = CON$Peak[CON_peak], 
                                    RetTime_CON = CON$RetTime[CON_peak],
                                    PeakArea_CON = CON$Area[CON_peak], 
                                    UV_Count_CON = FinalCount, 
                                    Subtracted_UV_Mean_CON = UV_Mean, 
                                    PeakRatio_CON = ratio))
        }   else  { 
            
            #The peak is declared NOT to have a match in the control and given NA
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[CC_peak], 
                                    RetTime_CC = Coculture$RetTime[CC_peak], 
                                    PeakArea_CC = Coculture$Area[CC_peak], 
                                    PeakNo_CON = NA, RetTime_CON = NA, 
                                    PeakArea_CON = NA, UV_Count_CON = NA, 
                                    Subtracted_UV_Mean_CON = NA, 
                                    PeakRatio_CON = NA))
        }   
        CC_peak <- CC_peak + 1
    }
    return(Coculture_df)
}

#############################################
#6.CON Consolidator: Removes double peak matching to a coculture peak
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
        }   else if (!is.na(Coculture_df[RowNo, 4]) && 
                     !is.na(Coculture_df[RowNo, 10] && 
                            abs(Coculture_df[RowNo, 8]) < 
                            abs(Coculture_df[RowNo, 14]))) { 
            
            #When there are two peaks matched (!is.na for both),
            #the UV_means are compared, with the higher removed.
            
            MatchedPeak_df[RowNo, 1] <- CON1_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
        }   else if (!is.na(Coculture_df[RowNo, 4]) && 
                     !is.na(Coculture_df[RowNo, 10] && 
                            abs(Coculture_df[RowNo, 8]) > 
                            abs(Coculture_df[RowNo, 14]))) { 
            MatchedPeak_df[RowNo, 1] <- CON2_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
        }   else if (is.na(Coculture_df[RowNo, 10]))   {
            
            #If no double-peak mactching but signifies a matched CON1 peak
            #then CON1 peak set as the matched peak
            
            MatchedPeak_df[RowNo, 1] <- CON1_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
        }   else {
            
            #If no double-peak mactching or matched CON1 peak
            #then CON2 peak set as the matched peak
            
            MatchedPeak_df[RowNo, 1] <- CON2_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
        }
        RowNo <- RowNo + 1
    }
    return(MatchedPeak_df)
}

#############################################
#7.Metabolite Effect Characteriser: Converting PeakRatio to a Factor
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
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] > -100
                     && Refined_Coculture_df[RowNo, 11] <= -20) {
            Metabolite_effect_df[RowNo, 1] <- 2 #Suppression
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] > -20
                     && Refined_Coculture_df[RowNo, 11] < 20) {
            Metabolite_effect_df[RowNo, 1] <- 3 #Little to No Change
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] >= 20
                     && Refined_Coculture_df[RowNo, 11] < 100) {
            Metabolite_effect_df[RowNo, 1] <- 4 #Enhancement
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] >= 100) {
            Metabolite_effect_df[RowNo, 1] <- 5 #Major Enhancement
        }
        RowNo <- RowNo + 1
    }
    Refined_Coculture_df <- 
        cbind(Refined_Coculture_df, Metabolite_effect_df)
    return(Refined_Coculture_df)
}

#############################################
#MIMI_Part_1: The main working-function to compare the peak-matching and refinement
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
        
        #Runs the peak matching algorithm for CON1 and CON2 separately.
        #Then merges the two together into a unified df.
        
        Coculture_df <- Peak_Matcher(CON1, Coculture, CON1_UV, Coculture_UV)
        Coculture_df2 <- Peak_Matcher(CON2, Coculture, CON2_UV, Coculture_UV)
        Coculture_df_merged <- cbind(Coculture_df, Coculture_df2[, 4:9])
        
        #Removes the first arbitrary row of missing values
        
        Coculture_df_merged <- Coculture_df_merged[2:nrow(Coculture_df_merged), ]
        
        #Performs a function to correct for double peak matching to a unique
        #peak to peaks from more than one CON
        
        MatchedPeak_df <- CON_Consolidator(Coculture_df_merged, CON1_Name, CON2_Name)
        
        #Final processing steps in creating the tidied df (Refined_Coculture_df)
        
        MatchedPeak_df <- MatchedPeak_df[1:nrow(MatchedPeak_df), ]
        df <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(Coculture_df_merged))), 
                       "Sample_Ref")
        df[1:nrow(df), ] <- Coculture_Name
        Refined_Coculture_df <- Coculture_df[1:nrow(Coculture_df_merged), 1:3]
        Refined_Coculture_df <- cbind(df, Refined_Coculture_df)
        Refined_Coculture_df <- cbind(Refined_Coculture_df, MatchedPeak_df)
        
        #Performs a function to characterise effects onto metabolites
        #and adds this into the tidied data set
        
        Refined_Coculture_df <- 
            Effect_Categoriser(Refined_Coculture_df, Coculture_Name)
        
        #Rewrites the tidied dataset to a csv file
        
        write.csv(Refined_Coculture_df, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         Coculture_Name, ".CSV"), row.names = FALSE)
    }
}

#Load in excel file named Interaction_Matrix and use MIMI() function:
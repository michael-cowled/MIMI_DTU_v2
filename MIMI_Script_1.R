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
#1.ReadExcel
#2.UVcheck
#3.ReadUV
#4.UVSubtract
#5.PeakMatcher
#6.CONConsolidator
#7.EffectCategoriser

#Special Terms:
#1. CC - Coculture: the duoculture of two organisms grown together.
#2. CON - Control: the axenic monoculture of one of the organisms.

#Reading in the functions:

#############################################
#1.ReadExcel: Creates 3 data frames for the 3 samples to be compared.
#############################################

ReadExcel <- function(Excel.Name) {
    
    #Input is Excel.Name which is read from the Interaction_Matrix object
    
    #Read in a df based on NovaC excel format
    raw.df <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", 
                                Excel.Name, ".xlsm"), skip = 3)
    
    #A new df is set up to capture the top 5 UV maxima for each peak
    UV.df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
                      c("UV1", "UV2", "UV3", "UV4", "UV5"))
    
    #The following extracts and sorts the top 5 UV maxima from each peak
    n <- nrow(raw.df)
    CC.peak <- 1
    while (CC.peak <= n) {
        if (!is.na(raw.df$'UV Peaks')) {
            raw.UV <-raw.df$`UV Peaks`[[CC.peak]]
            raw.UV.processed <- gsub("\\(","", raw.UV)
            raw.UV.processed <- gsub("\\)","", raw.UV.processed)
            raw.UV.processed <- gsub("< 190","", raw.UV.processed)
            raw.UV.processed <- gsub("s","", raw.UV.processed)
            raw.UV.split <- strsplit(raw.UV.processed, "\r\n")
            raw.UV.split <- as.data.frame(raw.UV.split, col.names = "UV1")
            raw.UV.split <- separate(raw.UV.split, 'UV1', c("UV1", "UV1_percent"), 
                                     sep = " ")
            raw.UV.split <- transform(raw.UV.split, UV1_percent = as.numeric(UV1_percent))
            raw.UV.ordered <- raw.UV.split[order(raw.UV.split$UV1_percent, decreasing = TRUE),]
            
            #A second UV.df is set up to capture the UV data for a single peak 
            UV.df.peak <- setNames(data.frame(matrix(ncol = 5, nrow = 1)), 
                                   c("UV1", "UV2", "UV3", "UV4", "UV5"))
            
            #Transposes UVs into row format; more efficient method possible
            UV.df.peak[1,1] <- raw.UV.ordered[1,1]
            UV.df.peak[1,2] <- raw.UV.ordered[2,1]
            UV.df.peak[1,3] <- raw.UV.ordered[3,1]
            UV.df.peak[1,4] <- raw.UV.ordered[4,1]
            UV.df.peak[1,5] <- raw.UV.ordered[5,1]
            
            #Lists the UV for the associated peak into the originally setup UV.df
            UV.df <- rbind(UV.df, UV.df.peak)
            
        }   else {
            
            #Lists the UVs as nothing instead
            UV.df.peak <- c(NA, NA, NA, NA, NA)
            UV.df <- rbind(UV.df, UV.df.peak)
            
        }   
        CC.peak <- CC.peak + 1
    }
    
    #UV.df is added onto the loaded raw.df object, correlating peaks with UV
    names(UV.df) <- c("UV1", "UV2", "UV3", "UV4", "UV5")   
    UV.df <- transform(UV.df, UV1 = as.numeric(UV1))
    UV.df <- transform(UV.df, UV2 = as.numeric(UV2))
    UV.df <- transform(UV.df, UV3 = as.numeric(UV3))
    UV.df <- transform(UV.df, UV4 = as.numeric(UV4))
    UV.df <- transform(UV.df, UV5 = as.numeric(UV5))
    raw.df.with.sorted.UV <- cbind(raw.df, UV.df)
    raw.df.with.sorted.UV <- select(raw.df.with.sorted.UV, 1:4, 20:24)
    return(raw.df.with.sorted.UV)
    
}

#############################################
#2.UV_Check: Verifies the matching UV maxima for peaks in CON1/2 & Coculture
#############################################

#UVcheck: Comparing CON1 to coculture
#Inputs are CON1/CON2 and coculture as read from the Interaction_Matrix
#'CC.peak' refers to the peak no. being compared in Coculture
#'CON.peak' refers to the peak no. being compared in control
#UV data is located in columns 5:9, hence 'CON.UV.no' and 'CC.UV.no' set to start at 5.

UVCheck <- function(control, Coculture, CC.peak, CON.peak) {
    
    UV.count <- 0
    CC.UV.no <- 5 
    CON.UV.no <- 5
    
    while (CC.UV.no < 10) {
        CON.UV.no <- 5
        while (CON.UV.no < 10 && CC.UV.no < 10) {
            if (CON.UV.no == 9)  {
                CON.UV.no <- CON.UV.no + 1
                CC.UV.no <- CC.UV.no + 1
            }   else if (is.na(Coculture[CC.peak,CC.UV.no])) {
                CC.UV.no <- CC.UV.no + 1
            }   else if (is.na(control[CON.peak,CON.UV.no])) {
                CON.UV.no <- CON.UV.no + 1
            }   else if (Coculture[CC.peak,CC.UV.no] <= control[CON.peak,CON.UV.no] + 2 && 
                         Coculture[CC.peak,CC.UV.no] >= control[CON.peak,CON.UV.no] - 2) {
                UV.count <- UV.count + 1
                CON.UV.no <- CON.UV.no + 1
            }   else {
                CON.UV.no <- CON.UV.no + 1
            }
        }
    }
    return(UV.count)
}

#############################################
#3.ReadUV Function: Reads in the raw UV data for every wavelength
#############################################

#Input is Excel.Name which is read from the Interaction_Matrix object
#Generates a df with col1 = Wavelength, col2 = Abs for Peak 1, etc.

ReadUV <- function(Excel.Name) {
    UV.Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/",
                                 Excel.Name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(UV.Name)
}

#############################################
#4.UVSubtract: Subtracts whole UV spectra and calculates the mean for peaks 
#in CON1/2 & Coculture
#############################################

#UVsubtract: Comparing control to coculture

UVSubtract <- function(CON.UV, CC.UV, CC.peak, CON.peak) {
    
    #Inputs are CON1.UV and CC.UV as created from the Interaction_Matrix
    #'CC.peak' refers to the peak no. being compared in Coculture
    #'CON.peak' refers to the peak no. being compared in CONtrol
    
    CC.UV.no <- CC.peak + 1
    CON.UV.no <- CON.peak + 1
    
    #Subtracts absorbances at each wavelength and calculates the column mean
    #'CC.UV.no' refers to the column of abs data being compared in Coculture
    #'CON.UV.no' refers to the column of abs data being compared in control
    #Note: Col2 = Peak#1, Col3 = Peak#2, etc. hence '+1'
    
    SubtractedUV <- CC.UV[,CC.UV.no]-CON.UV[,CON.UV.no]
    UV.Mean <- mean(SubtractedUV)
    return(UV.Mean)
}

#############################################
#5.PeakMatcher: Finds and compares the nearest matching peak in control
#############################################

PeakMatcher <- function(CON, Coculture, CON.UV, CC.UV) {
    
    #Sets up a df to with the desired column names.
    CC.df <- setNames(data.frame(matrix(ncol = 9, nrow = 0)), 
                             c("PeakNo_CC", "RetTime_CC", "PeakArea_CC", "PeakNo_CON", 
                               "RetTime_CON", "PeakArea_CON", "UV_Count_CON",
                               "Subtracted_UV_Mean_CON", "PeakRatio_CON"))
    
    n <- nrow(Coculture)
    CC.peak <- 1
    
    #'CC.peak' corresponds to the peak no. to be compared in the coculture
    
    while (CC.peak < n + 1) {
        CON.peak <- which(abs(CON$RetTime-Coculture$RetTime[CC.peak]) ==
                              min(abs(CON$RetTime-Coculture$RetTime[CC.peak])))
        
        #'CON.peak' finds the peak in CON1 with the closest ret-time to peak 'CC.peak'
        
        ratio = (((Coculture[CC.peak,3] - CON[CON.peak,3])/CON[CON.peak,3])*100) 
        
        #Computes the ratio of peak areas as a %
        
        Final.Count <- UVcheck(CON, Coculture, CC.peak, CON.peak)
        UV.Mean <- UVSubtract(CON.UV, CC.UV, CC.peak, CON.peak)
        
        #Performs both the UVcheck and UVsubtract functions
        
        if (Coculture$RetTime[CC.peak] < (CON$RetTime[CON.peak] + 0.15) && 
            Coculture$RetTime[CC.peak] > (CON$RetTime[CON.peak] -0.15) &&
            (Final.Count > 2 || abs(UV.Mean) < 1.5)) {
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UV.count is at least 3 OR UV.Mean < 1.5
            #If satisfied, the peak is declared a match and added to df
            
            CC.df <- rbind(CC.df, 
                                  c(PeakNo_CC = Coculture$Peak[CC.peak], 
                                    RetTime_CC = Coculture$RetTime[CC.peak], 
                                    PeakArea_CC = Coculture$Area[CC.peak], 
                                    PeakNo_CON = CON$Peak[CON.peak], 
                                    RetTime_CON = CON$RetTime[CON.peak],
                                    PeakArea_CON = CON$Area[CON.peak], 
                                    UV_Count_CON = Final.Count, 
                                    Subtracted_UV_Mean_CON = UV.Mean, 
                                    PeakRatio_CON = ratio))
            
        }    else if (Coculture$RetTime[CC.peak] < (CON$RetTime[CON.peak] + 0.15) && 
                      Coculture$RetTime[CC.peak] > (CON$RetTime[CON.peak] -0.15) &&
                      Final.Count > 1 && abs(UV.Mean) < 2) {  
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UV.count is at least 2 AND UV.Mean < 2
            #Both UV.count and UV.Mean are less strict requirement
            #But requiring satisfaction of both adds stringency
            #If satisfied, the peak is declared a match and added to df
            
            CC.df <- rbind(CC.df, 
                                  c(PeakNo_CC = Coculture$Peak[CC.peak], 
                                    RetTime_CC = Coculture$RetTime[CC.peak], 
                                    PeakArea_CC = Coculture$Area[CC.peak], 
                                    PeakNo_CON = CON$Peak[CON.peak], 
                                    RetTime_CON = CON$RetTime[CON.peak],
                                    PeakArea_CON = CON$Area[CON.peak], 
                                    UV_Count_CON = Final.Count, 
                                    Subtracted_UV_Mean_CON = UV.Mean, 
                                    PeakRatio_CON = ratio))
            
        }    else if (Coculture$RetTime[CC.peak] < (CON$RetTime[CON.peak] + 0.05) && 
                      Coculture$RetTime[CC.peak] > (CON$RetTime[CON.peak] -0.05) &&
                      Final.Count > 1) {  
            
            CC.df <- rbind(CC.df, 
                                  c(PeakNo_CC = Coculture$Peak[CC.peak], 
                                    RetTime_CC = Coculture$RetTime[CC.peak], 
                                    PeakArea_CC = Coculture$Area[CC.peak], 
                                    PeakNo_CON = CON$Peak[CON.peak], 
                                    RetTime_CON = CON$RetTime[CON.peak],
                                    PeakArea_CON = CON$Area[CON.peak], 
                                    UV_Count_CON = Final.Count, 
                                    Subtracted_UV_Mean_CON = UV.Mean, 
                                    PeakRatio_CON = ratio))
        }   else  { 
            
            #The peak is declared NOT to have a match in the control and given NA
            
            CC.df <- rbind(CC.df, 
                                  c(PeakNo_CC = Coculture$Peak[CC.peak], 
                                    RetTime_CC = Coculture$RetTime[CC.peak], 
                                    PeakArea_CC = Coculture$Area[CC.peak], 
                                    PeakNo_CON = NA, RetTime_CON = NA, 
                                    PeakArea_CON = NA, UV_Count_CON = NA, 
                                    Subtracted_UV_Mean_CON = NA, 
                                    PeakRatio_CON = NA))
        }   
        CC.peak <- CC.peak + 1
    }
    return(CC.df)
}

#############################################
#6.CONConsolidator: Removes double peak matching to a coculture peak
#############################################

#The following code removes double assignments to a coculture peak
#from two different controls.

#An example is that perhaps peak 1 from the coculture matches
#to peak 1 in CON1 and peak 2 in CON2.

#This function will determine which control's matched peak matches
#closest and remove the assignment for the weakest match.

CONConsolidator <- function(CC.df, CON1.Name, CON2.Name) {
    
    Row.no <- 1
    Total.Rows <- nrow(CC.df)
    
    #sets up a new tidier version of the df
    
    Matched.Peak.df <- setNames(data.frame(matrix(ncol = 7, 
                                                 nrow = nrow(CC.df))), 
                               c("Matched_CON", "PeakNo_CON", "RetTime_CON", 
                                 "PeakArea_CON", "UV_Count", 
                                 "Subtracted_UV.Mean", "PeakRatio"))
    
    while (Row.no <= Total.Rows + 1) {
        if (is.na(CC.df[Row.no, 4]) && 
            is.na(CC.df[Row.no, 10])) {
        }   else if (!is.na(CC.df[Row.no, 4]) && 
                     !is.na(CC.df[Row.no, 10] && 
                            abs(CC.df[Row.no, 8]) < 
                            abs(CC.df[Row.no, 14]))) { 
            
            #When there are two peaks matched (!is.na for both),
            #the UV.Means are compared, with the higher removed.
            
            Matched.Peak.df[Row.no, 1] <- CON1.Name
            Matched.Peak.df[Row.no, 2:7] <- CC.df[Row.no, 4:9]
        }   else if (!is.na(CC.df[Row.no, 4]) && 
                     !is.na(CC.df[Row.no, 10] && 
                            abs(CC.df[Row.no, 8]) > 
                            abs(CC.df[Row.no, 14]))) { 
            Matched.Peak.df[Row.no, 1] <- CON2.Name
            Matched.Peak.df[Row.no, 2:7] <- CC.df[Row.no, 10:15]
        }   else if (is.na(CC.df[Row.no, 10]))   {
            
            #If no double-peak mactching but signifies a matched CON1 peak
            #then CON1 peak set as the matched peak
            
            Matched.Peak.df[Row.no, 1] <- CON1.Name
            Matched.Peak.df[Row.no, 2:7] <- CC.df[Row.no, 4:9]
        }   else {
            
            #If no double-peak mactching or matched CON1 peak
            #then CON2 peak set as the matched peak
            
            Matched.Peak.df[Row.no, 1] <- CON2.Name
            Matched.Peak.df[Row.no, 2:7] <- CC.df[Row.no, 10:15]
        }
        Row.no <- Row.no + 1
    }
    return(Matched.Peak.df)
}

#############################################
#7.Metabolite Effect Characteriser: Converting PeakRatio to a Factor
#############################################

#The following piece of code adds a column that categorises the peak 
#areas into suppressions, and enhancements

EffectCategoriser <- function(Refined.CC.df, CC.Name) {
    
    Row.no <- 1
    Total.Rows <- nrow(Refined.CC.df)
    
    #A df is created to list the effects corresponding to matched peaks
    
    Metabolite.Effect.df <- setNames(data.frame(matrix(ncol=1, nrow=Total.Rows)), 
                                     c("Metabolite_Effect"))
    
    while (Row.no <= Total.Rows) {
        if (Refined.CC.df[Row.no, 1] == CC.Name && 
            is.na(Refined.CC.df[Row.no, 11])) {
            Metabolite.Effect.df[Row.no, 1] <- 6 #Induction
        }   else if (Refined.CC.df[Row.no, 1] == CC.Name && 
                     Refined.CC.df[Row.no, 11] > -100
                     && Refined.CC.df[Row.no, 11] <= -20) {
            Metabolite.Effect.df[Row.no, 1] <- 2 #Suppression
        }   else if (Refined.CC.df[Row.no, 1] == CC.Name && 
                     Refined.CC.df[Row.no, 11] > -20
                     && Refined.CC.df[Row.no, 11] < 20) {
            Metabolite.Effect.df[Row.no, 1] <- 3 #Little to No Change
        }   else if (Refined.CC.df[Row.no, 1] == CC.Name && 
                     Refined.CC.df[Row.no, 11] >= 20
                     && Refined.CC.df[Row.no, 11] < 100) {
            Metabolite.Effect.df[Row.no, 1] <- 4 #Enhancement
        }   else if (Refined.CC.df[Row.no, 1] == CC.Name && 
                     Refined.CC.df[Row.no, 11] >= 100) {
            Metabolite.Effect.df[Row.no, 1] <- 5 #Major Enhancement
        }
        Row.no <- Row.no + 1
    }
    Refined.CC.df <- 
        cbind(Refined.CC.df, Metabolite.Effect.df)
    return(Refined.CC.df)
}

#############################################
#MIMI_Part_1: The main working-function to compare the peak-matching and refinement
#############################################

MIMI <- function() {           
    
    #Will not work unless Interaction_Matrix object created by the user.
    
    Matrix.Total.Rows <- nrow(Interaction_Matrix)
    Matrix.Row.no <- 1
    
    while (Matrix.Row.no <= Matrix.Total.Rows) {
        
        CON1.Name <- as.character(Interaction_Matrix[Matrix.Row.no,1])
        CON2.Name <- as.character(Interaction_Matrix[Matrix.Row.no,2])
        CC.Name <- as.character(Interaction_Matrix[Matrix.Row.no,3])
        Matrix.Row.no <- Matrix.Row.no +1
        
        #Generates 3 dataframes using the ReadExcel function for the first
        #interction to be investigated from the Interaction_Matrix
        
        CON1 <- as.data.frame(ReadExcel(CON1.Name))
        CON2 <- as.data.frame(ReadExcel(CON2.Name))
        Coculture <- as.data.frame(ReadExcel(CC.Name))
        
        #Generates 3 dataframes using the ReadUV function
        
        CON1.UV <- as.data.frame(ReadUV(CON1.Name))
        CON2.UV <- as.data.frame(ReadUV(CON2.Name))
        CC.UV <- as.data.frame(ReadUV(CC.Name))
        
        #Runs the peak matching algorithm for CON1 and CON2 separately.
        #Then merges the two together into a unified df.
        
        CC.df <- PeakMatcher(CON1, Coculture, CON1.UV, CC.UV)
        CC.df2 <- PeakMatcher(CON2, Coculture, CON2.UV, CC.UV)
        CC.df.merged <- cbind(CC.df, CC.df2[, 4:9])
        
        #Performs a function to correct for double peak matching to a unique
        #peak to peaks from more than one CON
        
        Matched.Peak.df <- CONConsolidator(CC.df.merged, CON1.Name, CON2.Name)
        
        #Final processing steps in creating the tidied df (Refined.CC.df)
        
        Matched.Peak.df <- Matched.Peak.df[1:nrow(Matched.Peak.df), ]
        df <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(CC.df.merged))), 
                       "Sample_Ref")
        df[1:nrow(df), ] <- CC.Name
        Refined.CC.df <- CC.df[1:nrow(CC.df.merged), 1:3]
        Refined.CC.df <- cbind(df, Refined.CC.df)
        Refined.CC.df <- cbind(Refined.CC.df, Matched.Peak.df)
        
        #Performs a function to characterise effects onto metabolites
        #and adds this into the tidied data set
        
        Refined.CC.df <- 
            EffectCategoriser(Refined.CC.df, CC.Name)
        names(Refined.CC.df)[2:4] <- c("PeakNo_CC", "RetTime_CC", "PeakArea_CC")
        
        #Rewrites the tidied dataset to a csv file
        
        write.csv(Refined.CC.df, 
                  paste0("Testing Broad-Scale Interactions/OutputFiles/", 
                         CC.Name, ".CSV"), row.names = FALSE)
    }
}

#Load in excel file named Interaction_Matrix and use MIMI() function:
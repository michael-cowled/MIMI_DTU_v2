#Functions to be pre-loaded prior to MIMI()

library(dplyr)
library(tidyr)
library(readxl)

#Double-Peak Remover
Double_Peak_Remover <- function(Interaction_Matrix) {
    
    #To clean up data with multiple matches for a peak assignment in a control.
    
    #Generation of Logic_Table
    #To improve efficiency, a quick check is conducted first and a new TRUE/FALSE table produced in order to filter.
    
    Logic_TotalRows <- 1
    
    while (Logic_TotalRows > 0) {
        
        Matrix_TotalRows <- nrow(Interaction_Matrix)
        Matrix_Row_No <- 1
        Logic_Table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("Coculture_Name"))
        
        while (Matrix_Row_No <= Matrix_TotalRows) {
            
            Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
            Coculture_Name2 <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), c("Coculture_Name"))
            Coculture_Name2[1,1] <- Coculture_Name
            df_Name <- read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".csv"))
            df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), sep = "-", remove = FALSE)
            df_Name <- filter(df_Name, Combined !="NA-NA")
            logic <- length(unique(df_Name$Combined)) == nrow(df_Name)
            if (logic == TRUE) {
                Matrix_Row_No <- Matrix_Row_No +1
            } else {
                Logic_Table <- rbind(Logic_Table, Coculture_Name2)
                Matrix_Row_No <- Matrix_Row_No +1
            }
        }
        
        #Now to go through and fix the values in the logic table, perhaps have a final check to go back through in a loop
        
        Logic_TotalRows <- nrow(Logic_Table)
        Logic_Row_No <- 1
        
        while (Logic_Row_No <= Logic_TotalRows) {
            Coculture_Name <- as.character(Logic_Table[Logic_Row_No,1])
            print(Coculture_Name)
            df_Name <- read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".csv"))
            df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), sep = "-", remove = FALSE)
            df_Name4 <- df_Name
            df_Name <- filter(df_Name, Combined != "NA-NA")
            df_Name$Duplicated <- duplicated(df_Name$Combined)
            df_Name2 <- filter(df_Name, Duplicated == TRUE)
            df_Name3 <- filter(df_Name, Combined == df_Name2[1, 5])
            if (df_Name3[1,10] > df_Name3[2,10]) {
                Bad_Peak <- df_Name3[2,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }    else if (df_Name3[1,10] < df_Name3[2,10]) {
                Bad_Peak <- df_Name3[1,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }   else if (df_Name3[1,11] < df_Name3[2,11]) {
                Bad_Peak <- df_Name3[2,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }   else {
                Bad_Peak <- df_Name3[1,2]
                df_Name4[Bad_Peak, 5:12] <- NA
            }
            df_Name4 <- select(df_Name4, -Combined)
            write.csv(df_Name4, paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".CSV"), row.names = FALSE)
            Logic_Row_No <- Logic_Row_No + 1
        }
        print("loop finished")    
    }
}

#Missing Control Peaks:

Missing_Control_Peaks <- function(Interaction_Matrix) {
    
    #Below code adds in the non-matched unique peaks from the controls.
    
    Matrix_TotalRows <- nrow(Interaction_Matrix)
    Matrix_Row_No <- 1
    
    while (Matrix_Row_No <= Matrix_TotalRows) {
        
        CON1_Name <- as.character(Interaction_Matrix[Matrix_Row_No,1])
        CON2_Name <- as.character(Interaction_Matrix[Matrix_Row_No,2])
        Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
        df_Name <- read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".csv"))
        df_Name$Sample_Ref <- as.character(df_Name$Sample_Ref)
        df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), sep = "-", remove = FALSE)
        CON1 <- as.data.frame(Read_Excel(CON1_Name))
        CON2 <- as.data.frame(Read_Excel(CON2_Name))
        
        Matrix_Row_No <- Matrix_Row_No + 1
        n <- nrow(CON1)
        i <- 1
        while (i <= n) {
            cat(CON1_Name, "-", i, "\n")
            if (any(df_Name[,5] == paste0(CON1_Name, "-", i), na.rm = TRUE)) {
                i <- i +1    
            }   else if (any(df_Name[,5] != paste0(CON1_Name, "-", i), na.rm = TRUE)) {
                df_Name <- rbind(df_Name, c(Sample_Ref = CON1_Name, PeakNo_CC = NA, RetTime_CC = NA, PeakArea_CC = NA,
                                            Combined = NA, MAtched_CON = NA, PeakNo_CON = CON1$Peak[i], 
                                            RetTime_CON = CON1$RetTime[i], PeakArea_CON = CON1$Area[i], 
                                            UV_Count = NA, Subtracted_UV_Mean = NA, PeakRatio = -100,
                                            Metabolite_Effect = 1))
                i <- i + 1
            }   else {
                i <- i + 1
            }
        }
        
        Matrix_Row_No <- Matrix_Row_No + 1
        n <- nrow(CON2)
        i <- 1
        while (i <= n) {
            cat(CON2_Name, "-", i, "\n")
            if (any(df_Name[,5] == paste0(CON2_Name, "-", i), na.rm = TRUE)) {
                i <- i +1    
            }   else if (any(df_Name[,5] != paste0(CON2_Name, "-", i), na.rm = TRUE)) {
                df_Name <- rbind(df_Name, c(Sample_Ref = CON2_Name, PeakNo_CC = NA, RetTime_CC = NA, PeakArea_CC = NA,
                                            Combined = NA, MAtched_CON = NA, PeakNo_CON = CON2$Peak[i], 
                                            RetTime_CON = CON2$RetTime[i], PeakArea_CON = CON2$Area[i], 
                                            UV_Count = NA, Subtracted_UV_Mean = NA, PeakRatio = -100,
                                            Metabolite_Effect = 1))
                i <- i + 1
            }   else {
                i <- i + 1
            }
        }
        df_Name <- select(df_Name, -Combined)
        write.csv(df_Name, paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".CSV"), row.names = FALSE)
        print("loop finished")
    }
}

#Read_Excel Function to create 3 data frames for the 3 samples to be compared.
Read_Excel <- function(Excel_Name) {
    df_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", Excel_Name, ".xlsm"), 
                          skip = 3)
    
    #The following extracts and sorts the top 5 UV maxima from each peak
    UV_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("UV1", "UV2", "UV3", "UV4", "UV5"))
    
    n <- nrow(df_Name)
    i <- 1
    while (i <= n) {
        if (!is.na(df_Name$'UV Peaks')) {
            UV_separate <-df_Name$`UV Peaks`[[i]] ##Looking at row by row. starting with compound 1
            UV_separate <- gsub("\\(","", UV_separate)
            UV_separate <- gsub("\\)","", UV_separate)
            UV_separate <- gsub("< 190","", UV_separate)
            UV_separate <- gsub("s","", UV_separate)
            UV_separate <- strsplit(UV_separate, "\r\n")
            UV_separate <- as.data.frame(UV_separate, col.names = "UV1")
            UV_separate <- separate(UV_separate, 'UV1', c("UV1", "UV1_percent"), sep = " ")
            UV_separate <- UV_separate[order(UV_separate$UV1_percent, decreasing =TRUE),]
            UV_df <- rbind(UV_df, UV_separate[1:5,1])
            i <- i +1    
        }   else {
            i <- i + 1
            UV_separate <- c(NA, NA, NA, NA, NA)
            UV_df <- rbind(UV_df, UV_separate[1:5,1])
        }
    }
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

#Read_UV Function
Read_UV <- function(Excel_Name) {
    UV_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", Excel_Name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(UV_Name)
}

#UV Check

#Test UVCheck Function for MIMI

UVcheck <- function(CON1, Coculture, i, z) {
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
            }   else if (is.na(CON1[z,j])) {
                j <- j+1
            }   else if (Coculture[i,k] <= CON1[z,j] + 2 && Coculture[i,k] >= CON1[z,j] - 2) {
                uvcount <- uvcount + 1
                j <- j+1
            }   else {
                j <- j+1
            }
        }
    }
    if (uvcount >= 1) {
        cat("metabolite confirmed with uv count =", uvcount, "with a ratio of", round(ratio, digits =0), "% \n")
        return(uvcount)
    }   else    {
        cat("but UV does not match \n")
        return(uvcount)
    }   
}

#For CON2

UVCheck2 <- function(CON2, Coculture, i, z) {
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
            }   else if (Coculture[i,k] <= CON2[z,j] + 2 && Coculture[i,k] >= CON2[z,j] - 2) {
                uvcount <- uvcount + 1
                j <- j+1
            }   else {
                j <- j+1
            }
        }
    }
    if (uvcount >= 1) {
        cat("metabolite confirmed with uv count =", uvcount, "with a ratio of", round(ratio, digits =0), "% \n")
        return(uvcount)
    }   else    {
        cat("but UV does not match \n")
        return(uvcount)
    }   
}

#UV Subtract

UVSubtract1 <- function(CON1_UV, Coculture_UV, i, z) {
    k <- i + 1
    j <- z + 1
    SubtractedUV <- Coculture_UV[,k]-CON1_UV[,j]
    UV_Mean <- mean(SubtractedUV)
    return(UV_Mean)
}

UVSubtract2 <- function(CON2_UV, Coculture_UV, i, z) {
    k <- i + 1
    j <- z + 1
    SubtractedUV <- Coculture_UV[,k]-CON2_UV[,j]
    UV_Mean <- mean(SubtractedUV)
    return(UV_Mean)
}

#MIMI- The main workable function, and the function to CALL.

MIMI <- function() {           
    #Name is Microbial Interaction Metabolite Integrator
    #Requires the user now to make an object called Interaction_Matrix as a df of CON1, CON2, and Coculture
    
    Matrix_TotalRows <- nrow(Interaction_Matrix)
    Matrix_Row_No <- 1
    while (Matrix_Row_No <= Matrix_TotalRows) {
        CON1_Name <- as.character(Interaction_Matrix[Matrix_Row_No,1])
        CON2_Name <- as.character(Interaction_Matrix[Matrix_Row_No,2])
        Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
        
        Matrix_Row_No <- Matrix_Row_No +1
        
        #Generates 6 dataframes based on the imput names in quotes
        CON1 <- as.data.frame(Read_Excel(CON1_Name))
        CON2 <- as.data.frame(Read_Excel(CON2_Name))
        Coculture <- as.data.frame(Read_Excel(Coculture_Name))
        CON1_UV <- as.data.frame(Read_UV(CON1_Name))
        CON2_UV <- as.data.frame(Read_UV(CON2_Name))
        Coculture_UV <- as.data.frame(Read_UV(Coculture_Name))
        
        Coculture_df <- data.frame(PeakNo_CC = NA, RetTime_CC = NA, PeakArea_CC = NA,
                                   PeakNo_CON1 = NA, RetTime_CON1 = NA, PeakArea_CON1 = NA,
                                   UV_Count_CON1 = NA, Subtracted_UV_Mean_CON1 = NA, PeakRatio_CON1 = NA)
        n <- nrow(Coculture)
        i = 1
        while (i < n+1) {
            z <- which(abs(CON1$RetTime-Coculture$RetTime[i])==min(abs(CON1$RetTime-Coculture$RetTime[i])))
            ratio = (((Coculture[i,3] - CON1[z,3])/CON1[z,3])*100) #Computes the ratio of peak areas as a %
            FinalCount <- UVcheck(CON1, Coculture, i, z)
            UV_Mean <- UVSubtract1(CON1_UV, Coculture_UV, i, z)
            if (Coculture$RetTime[i] < (CON1$RetTime[z] + 0.15) && Coculture$RetTime[i] > (CON1$RetTime[z] -0.15) 
                && (FinalCount > 2 || abs(UV_Mean) < 1.5)) {  
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON1$RetTime[z], digits =2), "min in the control ")
                Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                      PeakArea_CC = Coculture$Area[i], PeakNo_CON1 = CON1$Peak[z], 
                                                      RetTime_CON1 = CON1$RetTime[z], PeakArea_CON1 = CON1$Area[z], 
                                                      UV_Count_CON1 = FinalCount, Subtracted_UV_Mean_CON1 = UV_Mean, 
                                                      PeakRatio_CON1 = ratio))
                i <- i+1
            }    else if (Coculture$RetTime[i] < (CON1$RetTime[z] + 0.15) && Coculture$RetTime[i] > (CON1$RetTime[z] -0.15) 
                          && FinalCount > 1 && abs(UV_Mean) < 2) {  
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON1$RetTime[z], digits =2), "min in the control ")
                Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                      PeakArea_CC = Coculture$Area[i], PeakNo_CON1 = CON1$Peak[z], 
                                                      RetTime_CON1 = CON1$RetTime[z], PeakArea_CON1 = CON1$Area[z], 
                                                      UV_Count_CON1 = FinalCount, Subtracted_UV_Mean_CON1 = UV_Mean, 
                                                      PeakRatio_CON1 = ratio))
                i <- i+1
            }   else  { 
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min is NOT present in the control \n")
                Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                      PeakArea_CC = Coculture$Area[i], PeakNo_CON1 = NA, 
                                                      RetTime_CON1 = NA, PeakArea_CON1 = NA, 
                                                      UV_Count_CON1 = NA, Subtracted_UV_Mean_CON1 = NA, PeakRatio_CON1 = NA))
                i <- i+1
            }   
        }
        
        #Adding in section for testing against CON2
        
        Coculture_df2 <- data.frame(PeakNo_CC = NA, RetTime_CC = Coculture_Name, PeakArea_CC = NA,
                                    PeakNo_CON2 = NA, RetTime_CON2 = NA, PeakArea_CON2 = NA,
                                    UV_Count_CON2 = NA, Subtracted_UV_Mean_CON2 = NA, PeakRatio_CON2 = NA)
        n <- nrow(Coculture)
        i = 1
        while (i < n+1) {
            z <- which(abs(CON2$RetTime-Coculture$RetTime[i])==min(abs(CON2$RetTime-Coculture$RetTime[i])))
            ratio = (((Coculture[i,3] - CON2[z,3])/CON2[z,3])*100) #Computes the ratio of peak areas as a %
            FinalCount <- UVCheck2(CON2, Coculture, i, z)
            UV_Mean <- UVSubtract2(CON2_UV, Coculture_UV, i, z)
            if (Coculture$RetTime[i] < (CON2$RetTime[z] + 0.15) && Coculture$RetTime[i] > (CON2$RetTime[z] -0.15)
                && (FinalCount > 2 || abs(UV_Mean) < 1.5)) {  
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON2$RetTime[z], digits =2), "min in the control ")
                Coculture_df2 <- rbind(Coculture_df2, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                        PeakArea_CC = Coculture$Area[i], PeakNo_CON2 = CON2$Peak[z], 
                                                        RetTime_CON2 = CON2$RetTime[z], PeakArea_CON2 = CON2$Area[z], 
                                                        UV_Count_CON2 = FinalCount, Subtracted_UV_Mean_CON2 = UV_Mean, 
                                                        PeakRatio_CON2 = ratio))
                i <- i+1
            }    else if (Coculture$RetTime[i] < (CON2$RetTime[z] + 0.15) && Coculture$RetTime[i] > (CON2$RetTime[z] -0.15)
                          && FinalCount > 1 && abs(UV_Mean) < 2) {  
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON2$RetTime[z], digits =2), "min in the control ")
                Coculture_df2 <- rbind(Coculture_df2, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                        PeakArea_CC = Coculture$Area[i], PeakNo_CON2 = CON2$Peak[z], 
                                                        RetTime_CON2 = CON2$RetTime[z], PeakArea_CON2 = CON2$Area[z], 
                                                        UV_Count_CON2 = FinalCount, Subtracted_UV_Mean_CON2 = UV_Mean, 
                                                        PeakRatio_CON2 = ratio))
                i <- i+1
            }   else  { 
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min is NOT present in the control \n")
                Coculture_df2 <- rbind(Coculture_df2, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                        PeakArea_CC = Coculture$Area[i], PeakNo_CON2 = NA, 
                                                        RetTime_CON2 = NA, PeakArea_CON2 = NA, 
                                                        UV_Count_CON2 = NA, Subtracted_UV_Mean_CON2 = NA, PeakRatio_CON2 = NA))
                i <- i+1
            }   
        }
        Coculture_df <- cbind(Coculture_df, Coculture_df2[, 4:9])
        
        #The following piece of code converts Coculture_df into a dataframe matching only a single peak from each control.
        RowNo <- 1
        Rows <- nrow(Coculture_df)
        
        #The following code can be removed if and when the original MIMI code is changed to remove the NA row.
        Coculture_df <- Coculture_df[2:nrow(Coculture_df), ]
        
        MatchedPeak_df <- setNames(data.frame(matrix(ncol = 7, nrow = nrow(Coculture_df))), 
                                   c("Matched_CON", "PeakNo_CON", "RetTime_CON", "PeakArea_CON",
                                     "UV_Count", "Subtracted_UV_Mean", "PeakRatio"))
        
        while (RowNo <= Rows) {
            if (is.na(Coculture_df[RowNo, 4]) && is.na(Coculture_df[RowNo, 10])) {
                RowNo <- RowNo + 1
            }   else if (!is.na(Coculture_df[RowNo, 4]) && !is.na(Coculture_df[RowNo, 10] 
                                                                  && abs(Coculture_df[RowNo, 8]) < abs(Coculture_df[RowNo, 14]))) { #when two peaks are matched
                MatchedPeak_df[RowNo, 1] <- CON1_Name
                MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
                RowNo <- RowNo + 1
            }   else if (!is.na(Coculture_df[RowNo, 4]) && !is.na(Coculture_df[RowNo, 10] 
                                                                  && abs(Coculture_df[RowNo, 8]) > abs(Coculture_df[RowNo, 14]))) { 
                MatchedPeak_df[RowNo, 1] <- CON2_Name
                MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
                RowNo <- RowNo + 1
            }   else if (is.na(Coculture_df[RowNo, 10]))   {
                MatchedPeak_df[RowNo, 1] <- CON1_Name
                MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
                RowNo <- RowNo + 1
            }   else {
                MatchedPeak_df[RowNo, 1] <- CON2_Name
                MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
                RowNo <- RowNo + 1
            }
        }
        
        MatchedPeak_df <- MatchedPeak_df[1:nrow(MatchedPeak_df), ]
        df <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(Coculture_df))), "Sample_Ref")
        df[1:nrow(df), ] <- Coculture_Name
        Refined_Coculture_df <- Coculture_df[1:nrow(Coculture_df), 1:3] #The df containing Coculture data to be cbinded onto.
        Refined_Coculture_df <- cbind(df, Refined_Coculture_df)
        Refined_Coculture_df <- cbind(Refined_Coculture_df, MatchedPeak_df)
        
        #The following piece of code adds a column that categorises the peak areas into suppressions, and enhancements
        RowNo <- 1
        Rows <- nrow(Refined_Coculture_df)
        
        Metabolite_effect_df <- setNames(data.frame(matrix(ncol = 1, nrow = Rows)), c("Metabolite_Effect"))
        
        while (RowNo <= Rows) {
            if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && is.na(Refined_Coculture_df[RowNo, 11])) {
                Metabolite_effect_df[RowNo, 1] <- 6 #Induction
                RowNo <- RowNo + 1
            }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && Refined_Coculture_df[RowNo, 11] > -100
                         && Refined_Coculture_df[RowNo, 11] <= -20) {
                Metabolite_effect_df[RowNo, 1] <- 2 #Suppression
                RowNo <- RowNo + 1
            }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && Refined_Coculture_df[RowNo, 11] > -20
                         && Refined_Coculture_df[RowNo, 11] < 20) {
                Metabolite_effect_df[RowNo, 1] <- 3 #Little to No Change
                RowNo <- RowNo + 1
            }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && Refined_Coculture_df[RowNo, 11] >= 20
                         && Refined_Coculture_df[RowNo, 11] < 100) {
                Metabolite_effect_df[RowNo, 1] <- 4 #Enhancement
                RowNo <- RowNo + 1
            }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && Refined_Coculture_df[RowNo, 11] >= 100) {
                Metabolite_effect_df[RowNo, 1] <- 5 #Complete Suppression
                RowNo <- RowNo + 1
            }
        }
        
        Refined_Coculture_df <- cbind(Refined_Coculture_df, Metabolite_effect_df)
        
        write.csv(Refined_Coculture_df, paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".CSV"), row.names = FALSE)
        print("loop finished")
    }
    Double_Peak_Remover(Interaction_Matrix)
    Missing_Control_Peaks(Interaction_Matrix)
}

#Load in excel file named Interaction_Matrix and use MIMI() function:
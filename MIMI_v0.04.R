#Now to embed the previous MIMI function into this While loop

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
            if (Coculture$RetTime[i] < CON1$RetTime[z] + 0.2 && Coculture$RetTime[i] > CON1$RetTime[z] -0.2 
                && FinalCount > 2 || abs(UV_Mean) < 1.5) {  
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON1$RetTime[z], digits =2), "min in the control ")
                Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                      PeakArea_CC = Coculture$Area[i], PeakNo_CON1 = CON1$Peak[z], 
                                                      RetTime_CON1 = CON1$RetTime[z], PeakArea_CON1 = CON1$Area[z], 
                                                      UV_Count_CON1 = FinalCount, Subtracted_UV_Mean_CON1 = UV_Mean, 
                                                      PeakRatio_CON1 = ratio))
                i <- i+1
            }    else if (Coculture$RetTime[i] < CON1$RetTime[z] + 0.2 && Coculture$RetTime[i] > CON1$RetTime[z] -0.2 
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
            if (Coculture$RetTime[i] < CON2$RetTime[z] + 0.2 && Coculture$RetTime[i] > CON2$RetTime[z] -0.2
                && FinalCount > 2 || abs(UV_Mean) < 1.5) {  
                cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON2$RetTime[z], digits =2), "min in the control ")
                Coculture_df2 <- rbind(Coculture_df2, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                        PeakArea_CC = Coculture$Area[i], PeakNo_CON2 = CON2$Peak[z], 
                                                        RetTime_CON2 = CON2$RetTime[z], PeakArea_CON2 = CON2$Area[z], 
                                                        UV_Count_CON2 = FinalCount, Subtracted_UV_Mean_CON2 = UV_Mean, 
                                                        PeakRatio_CON2 = ratio))
                i <- i+1
            }    else if (Coculture$RetTime[i] < CON2$RetTime[z] + 0.2 && Coculture$RetTime[i] > CON2$RetTime[z] -0.2
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


MIMI()
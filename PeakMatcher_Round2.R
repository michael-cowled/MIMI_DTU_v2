Peak_Matcher_Round2 <- function(Interaction_Matrix) {
    
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
                (any(df_Name[,7] == CON1$RetTime[z1], na.rm = TRUE) |
                (any(df_Name[,7] == CON2$RetTime[z2], na.rm = TRUE))) {
            
            #Checks for peak matching already, and skips to the next peak.
                
                i <- i+1
            }   else if (Coculture$RetTime[i] < (CON1$RetTime[z1] + 0.05) && 
                         Coculture$RetTime[i] > (CON1$RetTime[z1] -0.05) &&
                         FinalCount1 > 1)    {
                
                #Checks if the closest match in CON1 satisfies this test.
                
                z <- z1
                ratio = (((Coculture[i,3] - CON1[z,3])/CON1[z,3])*100)
                
                #Assignments of the matched peak
                
                df_Name$Matched_CON[i] <- CON1_Name
                df_Name$PeakNo_CON[i] <- CON1$Peak[z]
                df_Name$RetTime_CON[i] <- CON1$RetTime[z]
                df_Name$PeakArea_CON[i] <- CON1$Area[z]
                df_Name$UV_Count[i] <- FinalCount1
                df_Name$PeakRatio[i] <- ratio

                i <- i+1
                
            }   else if (Coculture$RetTime[i] < (CON2$RetTime[z2] + 0.05) && 
                         Coculture$RetTime[i] > (CON2$RetTime[z2] -0.05) &&
                         FinalCount2 > 1)    {
                z <- z2
                ratio = (((Coculture[i,3] - CON2[z,3])/CON2[z,3])*100)
                df_Name$Matched_CON[i] <- CON2_Name
                df_Name$PeakNo_CON[i] <- CON2$Peak[z]
                df_Name$RetTime_CON[i] <- CON2$RetTime[z]
                df_Name$PeakArea_CON[i] <- CON2$Area[z]
                df_Name$UV_Count[i] <- FinalCount2
                df_Name$PeakRatio[i] <- ratio
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
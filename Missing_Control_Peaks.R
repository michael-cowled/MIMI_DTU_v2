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
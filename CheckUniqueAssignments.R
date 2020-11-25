#To clean up data with multple matches for a peak assignment in a control.
#To improve efficiency, a quick check is conducted first and a new TRUE/FALSE table produced in order to filter.

Matrix_TotalRows <- nrow(Interaction_Matrix)
Matrix_Row_No <- 1

while (Matrix_Row_No <= Matrix_TotalRows) {

    CON1_Name <- as.character(Interaction_Matrix[Matrix_Row_No,1])
    CON2_Name <- as.character(Interaction_Matrix[Matrix_Row_No,2])
    Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
    df_Name <- read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".csv"))
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
                                                    UV_Count = NA, Subtracted_UV_Mean = NA, PeakRatio = -100))
            print(df_Name)
            i <- i + 1
        }   else {
            i <- i + 1
            }
    }
    
    Matrix_Row_No <- Matrix_Row_No + 1
    n <- nrow(CON1)
    i <- 1
    while (i <= n) {
        cat(CON2_Name, "-", i, "\n")
        if (any(df_Name[,5] == paste0(CON2_Name, "-", i), na.rm = TRUE)) {
            i <- i +1    
        }   else if (any(df_Name[,5] != paste0(CON2_Name, "-", i), na.rm = TRUE)) {
            df_Name <- rbind(df_Name, c(Sample_Ref = CON2_Name, PeakNo_CC = NA, RetTime_CC = NA, PeakArea_CC = NA,
                                        Combined = NA, MAtched_CON = NA, PeakNo_CON = CON2$Peak[i], 
                                        RetTime_CON = CON2$RetTime[i], PeakArea_CON = CON2$Area[i], 
                                        UV_Count = NA, Subtracted_UV_Mean = NA, PeakRatio = -100))
            print(df_Name)
            i <- i + 1
        }   else {
            i <- i + 1
        }
    }
    write.csv(df_Name, paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".CSV"), row.names = FALSE)
    print("loop finished")
}



#### Above works for one control.


Logic_Table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), "Coculture_Name")
df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), sep = "-", remove = FALSE)
df_Name <- filter(df_Name, Combined !="NA-NA")
logic <- length(unique(df_Name$Combined)) == nrow(df_Name)

    if (logic == TRUE) {
        Matrix_Row_No <- Matrix_Row_No +1
    } else if (logic == FALSE) {
        Logic_Table <- rbind(Logic_Table, Coculture_Name = Coculture_Name)
        Matrix_Row_No <- Matrix_Row_No +1
    }
}

#Now to go through and fix the values in the logic table, perhaps have a final check to go back through in a loop

Logic_TotalRows <- nrow(Logic_Table)
Logic_Row_No <- 1

while (Logic_Row_No <= Logic_TotalRows) {
    
    Coculture_Name <- as.character(Logic_Table[Logic_Row_No,1])
    df_Name <- read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".csv"))
    df_Name <- unite(df_Name, Combined, c(Matched_CON, PeakNo_CON), sep = "-", remove = FALSE)
    row_i <- 1
    
    while (i < nrow(Logic_Table) + 1) {
        
        if (df_Name[i, 5] == "NA-NA")
            row_i <- row_i + 1
            if else (df_Name[i, 5])
            if (df_Name[i, 5] == (df_Name[j, 5]) {
            Logic_Row_No <- Logic_Row_No +1
        } else if (logic == FALSE) {
            Logic_Table <- rbind(Logic_Table, Coculture_Name = Coculture_Name)
            Matrix_Row_No <- Matrix_Row_No +1
        }
    }
}


#Draft code, probably doesn't work, return to when brain can work.
#Code is meant to find, compare and remove non-unique peak matches.

j <- 1
k <- j+1 
while (k < nrow(Logic_Table) + 1) {
    k <- j+1 
    while (j < nrow(Logic_Table) + 1 && k < nrow(Logic_Table) + 1) {
        if (j==nrow(Logic_Table))  {
            j <- j+1
            k <- k+1
        }   else if (df_Name[j, 5] == "NA-NA") {
            k <- k+1
        }   else if (df_Name[j, 5] == df_Name[k, 5] && df_Name[j, 10] > df_Name[k, 10]) {
            df_Name[k, 6:12] <- NA
            j <- j+1
        }   else if (df_Name[j, 5] == df_Name[k, 5] && df_Name[j, 10] < df_Name[k, 10]) {
            df_Name[j, 6:12] <- NA
            j <- j+1
        }   else if (df_Name[j, 5] == df_Name[k, 5] && abs(df_Name[j, 11]) < abs(df_Name[k, 11]) {
            df_Name[k, 6:12] <- NA
            j <- j+1
        }   else if (df_Name[j, 5] == df_Name[k, 5] && abs(df_Name[j, 11]) > abs(df_Name[k, 11]) {
            df_Name[j, 6:12] <- NA
            j <- j+1
        }   else if (df_Name[j, 5] != df_Name[k, 5]) {
            df_Name[j, 6:12] <- NA
            j <- j+1
        }   
    }
}

Double_Peak_Remover <- function(Interaction_Matrix) {

#To clean up data with multiple matches for a peak assignment in a control.

#Generation of Logic_Table
#To improve efficiency, a quick check is conducted first and a new TRUE/FALSE table produced in order to filter.

Logic_TotalRows <- 1

while (Logic_TotalRows > 0) {
    
    Matrix_TotalRows <- nrow(Interaction_Matrix)
    Matrix_Row_No <- 1
    Logic_Table <- setNames(data.frame(matrix(ncol = 1, nrow = 0)), "Coculture_Name")
    
    while (Matrix_Row_No <= Matrix_TotalRows) {
        
        Coculture_Name <- as.character(Interaction_Matrix[Matrix_Row_No,3])
        df_Name <- read.csv(paste0("Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".csv"))
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
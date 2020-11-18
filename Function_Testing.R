MIMI <- function(CON1_Name, CON2_Name, Coculture_Name) {           
    #Name is Microbial Interaction Metabolite Integrator
    #Requires the user to make variables for CON1_Name, CON2_Name and Coculture_Name
    
    #Generates 3 dataframes based on the imput names in quotes
    CON1 <- as.data.frame(Read_Excel(CON1_Name))
    CON2 <- as.data.frame(Read_Excel(CON2_Name))
    Coculture <- as.data.frame(Read_Excel(Coculture_Name))
    
    Coculture_df <- data.frame(PeakNo_CC = Coculture_Name, RetTime_CC = NA, PeakArea_CC = NA,
                               PeakNo_CON = NA, RetTime_CON = CON1_Name, PeakArea_CON = NA,
                               UV_Count = NA, PeakRatio = NA)
    n <- nrow(Coculture)
    i = 1
    while (i < n+1) {
        z <- which(abs(CON1$RetTime-Coculture$RetTime[i])==min(abs(CON1$RetTime-Coculture$RetTime[i])))
        ratio = ((Coculture[i,3]/CON1[z,3])*100) #Computes the ratio of peak areas as a %
        FinalCount <- UVcheck(CON1, Coculture, i, z)
        if (Coculture$RetTime[i] < CON1$RetTime[z] + 0.2 && Coculture$RetTime[i] > CON1$RetTime[z] -0.2 && FinalCount > 0) {  
            cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON1$RetTime[z], digits =2), "min in the control ")
            Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                  PeakArea_CC = Coculture$Area[i], PeakNo_CON = CON1$Peak[z], 
                                                  RetTime_CON = CON1$RetTime[z], PeakArea_CON = CON1$Area[z], 
                                                  UV_Count = FinalCount, PeakRatio = ratio))
            i <- i+1
        }   else  { 
            cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min is NOT present in the control \n")
            Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                  PeakArea_CC = Coculture$Area[i], PeakNo_CON = NA, 
                                                  RetTime_CON = NA, PeakArea_CON = NA, 
                                                  UV_Count = NA, PeakRatio = NA))
            i <- i+1
        }   
    }
    
    #Adding in section for testing against CON2
    
    Coculture_df2 <- data.frame(PeakNo_CC = Coculture_Name, RetTime_CC = NA, PeakArea_CC = NA,
                                PeakNo_CON = NA, RetTime_CON = CON2_Name, PeakArea_CON = NA,
                                UV_Count = NA, PeakRatio = NA)
    n <- nrow(Coculture)
    i = 1
    while (i < n+1) {
        z <- which(abs(CON2$RetTime-Coculture$RetTime[i])==min(abs(CON2$RetTime-Coculture$RetTime[i])))
        ratio = ((Coculture[i,3]/CON1[z,3])*100) #Computes the ratio of peak areas as a %
        FinalCount <- UVCheck2(CON2, Coculture, i, z)
        if (Coculture$RetTime[i] < CON2$RetTime[z] + 0.2 && Coculture$RetTime[i] > CON2$RetTime[z] -0.2 && FinalCount > 0) {  
            cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON2$RetTime[z], digits =2), "min in the control ")
            Coculture_df2 <- rbind(Coculture_df2, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                    PeakArea_CC = Coculture$Area[i], PeakNo_CON = CON2$Peak[z], 
                                                    RetTime_CON = CON2$RetTime[z], PeakArea_CON = CON2$Area[z], 
                                                    UV_Count = FinalCount, PeakRatio = ratio))
            i <- i+1
        }   else  { 
            cat("Peak#", i, "@", round(Coculture$RetTime[i], digits =2), "min is NOT present in the control \n")
            Coculture_df2 <- rbind(Coculture_df2, c(PeakNo_CC = Coculture$Peak[i], RetTime_CC = Coculture$RetTime[i], 
                                                    PeakArea_CC = Coculture$Area[i], PeakNo_CON = NA, 
                                                    RetTime_CON = NA, PeakArea_CON = NA, 
                                                    UV_Count = NA, PeakRatio = NA))
            i <- i+1
        }   
    }
    Coculture_df <- cbind(Coculture_df, Coculture_df2[, 4:8])
    write.csv(Coculture_df, paste0("C:/Users/Radicinol/Desktop/R/Testing Broad-Scale Interactions/OutputFiles/", Coculture_Name, ".CSV"), row.names = FALSE)
    print("loop finished")
}


MIMI("27CON", "04CON", "27v04")
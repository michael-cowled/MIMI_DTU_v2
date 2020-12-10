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
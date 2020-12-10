Peak_Matcher <- function(CON1, Coculture, CON1_UV, Coculture_UV) {
    
    Coculture_df <- data.frame(PeakNo_CC = NA, RetTime_CC = NA, 
                               PeakArea_CC = NA, PeakNo_CON1 = NA, 
                               RetTime_CON1 = NA, PeakArea_CON1 = NA,
                               UV_Count_CON1 = NA, 
                               Subtracted_UV_Mean_CON1 = NA, 
                               PeakRatio_CON1 = NA)
    
    n <- nrow(Coculture)
    i = 1
    
    #'i' corresponds to the peak no. to be compared in the coculture
    
    while (i < n+1) {
        z <- which(abs(CON1$RetTime-Coculture$RetTime[i]) ==
                       min(abs(CON1$RetTime-Coculture$RetTime[i])))
        
        #'z' finds the peak in CON1 with the closest ret-time to peak 'i'
        
        ratio = (((Coculture[i,3] - CON1[z,3])/CON1[z,3])*100) 
        
        #Computes the ratio of peak areas as a %
        
        FinalCount <- UVcheck(CON1, Coculture, i, z)
        UV_Mean <- UVSubtract1(CON1_UV, Coculture_UV, i, z)
        
        #Performs both the UVcheck and UVsubtract functions
        
        if (Coculture$RetTime[i] < (CON1$RetTime[z] + 0.15) && 
            Coculture$RetTime[i] > (CON1$RetTime[z] -0.15) &&
            (FinalCount > 2 || abs(UV_Mean) < 1.5)) {
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UVcount is at least 3 OR UV_Mean < 1.5
            #If satisfied, the peak is declared a match and added to df
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[i], 
                                    RetTime_CC = Coculture$RetTime[i], 
                                    PeakArea_CC = Coculture$Area[i], 
                                    PeakNo_CON1 = CON1$Peak[z], 
                                    RetTime_CON1 = CON1$RetTime[z],
                                    PeakArea_CON1 = CON1$Area[z], 
                                    UV_Count_CON1 = FinalCount, 
                                    Subtracted_UV_Mean_CON1 = UV_Mean, 
                                    PeakRatio_CON1 = ratio))
            i <- i+1
        }    else if (Coculture$RetTime[i] < (CON1$RetTime[z] + 0.15) && 
                      Coculture$RetTime[i] > (CON1$RetTime[z] -0.15) &&
                      FinalCount > 1 && abs(UV_Mean) < 2) {  
            
            #Checks if the ret.time is within 0.15 min of each other
            #Checks the UVcount is at least 2 AND UV_Mean < 2
            #Both UVcount and UV_Mean are less strict requirement
            #But requiring satisfaction of both adds stringency
            #If satisfied, the peak is declared a match and added to df
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[i], 
                                    RetTime_CC = Coculture$RetTime[i], 
                                    PeakArea_CC = Coculture$Area[i], 
                                    PeakNo_CON1 = CON1$Peak[z], 
                                    RetTime_CON1 = CON1$RetTime[z],
                                    PeakArea_CON1 = CON1$Area[z], 
                                    UV_Count_CON1 = FinalCount, 
                                    Subtracted_UV_Mean_CON1 = UV_Mean, 
                                    PeakRatio_CON1 = ratio))
            i <- i+1
        }   else  { 
            
            #The peak is declared NOT to have a match in CON1, and NA added
            
            Coculture_df <- rbind(Coculture_df, 
                                  c(PeakNo_CC = Coculture$Peak[i], 
                                    RetTime_CC = Coculture$RetTime[i], 
                                    PeakArea_CC = Coculture$Area[i], 
                                    PeakNo_CON1 = NA, RetTime_CON1 = NA, 
                                    PeakArea_CON1 = NA, UV_Count_CON1 = NA, 
                                    Subtracted_UV_Mean_CON1 = NA, 
                                    PeakRatio_CON1 = NA))
            i <- i+1
        }   
    }
    return(Coculture_df)
}
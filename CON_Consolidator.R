#The following code removes double assignments to a coculture peak
#from two different controls.

#An example is that perhaps peak 1 from the coculture matches
#to peak 1 in CON1 and peak 2 in CON2.

#This function will determine which control's matched peak matches
#closest and remove the assignment for the weakest match.

CON_Consolidator <- function(Coculture_df, CON1_Name, CON2_Name) {

    RowNo <- 1
    Rows <- nrow(Coculture_df)
    
    #sets up a new tidier version of the df
    
    MatchedPeak_df <- setNames(data.frame(matrix(ncol = 7, 
                                                 nrow = nrow(Coculture_df))), 
                               c("Matched_CON", "PeakNo_CON", "RetTime_CON", 
                                 "PeakArea_CON", "UV_Count", 
                                 "Subtracted_UV_Mean", "PeakRatio"))
    
    while (RowNo <= Rows+1) {
        if (is.na(Coculture_df[RowNo, 4]) && 
            is.na(Coculture_df[RowNo, 10])) {
            RowNo <- RowNo + 1
        }   else if (!is.na(Coculture_df[RowNo, 4]) && 
                     !is.na(Coculture_df[RowNo, 10] && 
                            abs(Coculture_df[RowNo, 8]) < 
                            abs(Coculture_df[RowNo, 14]))) { 
            
            #When there are two peaks matched (!is.na for both),
            #the UV_means are compared, with the higher removed.
            
            MatchedPeak_df[RowNo, 1] <- CON1_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
            RowNo <- RowNo + 1
        }   else if (!is.na(Coculture_df[RowNo, 4]) && 
                     !is.na(Coculture_df[RowNo, 10] && 
                            abs(Coculture_df[RowNo, 8]) > 
                            abs(Coculture_df[RowNo, 14]))) { 
            MatchedPeak_df[RowNo, 1] <- CON2_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
            RowNo <- RowNo + 1
        }   else if (is.na(Coculture_df[RowNo, 10]))   {
            
            #If no double-peak mactching but signifies a matched CON1 peak
            #then CON1 peak set as the matched peak
            
            MatchedPeak_df[RowNo, 1] <- CON1_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
            RowNo <- RowNo + 1
        }   else {
            
            #If no double-peak mactching or matched CON1 peak
            #then CON2 peak set as the matched peak
            
            MatchedPeak_df[RowNo, 1] <- CON2_Name
            MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
            RowNo <- RowNo + 1
        }
    }
    return(MatchedPeak_df)
}
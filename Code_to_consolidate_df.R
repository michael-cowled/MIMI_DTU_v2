#remove this later
Coculture_df <- read.csv("Testing Broad-Scale Interactions/OutputFiles/27v04.CSV")

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
        print(RowNo)
    }   else if (!is.na(Coculture_df[RowNo, 4]) && !is.na(Coculture_df[RowNo, 10] 
                                                          && abs(Coculture_df[RowNo, 8]) < abs(Coculture_df[RowNo, 14]))) { #when two peaks are matched
        MatchedPeak_df[RowNo, 1] <- CON1_Name
        MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
        RowNo <- RowNo + 1
        print(RowNo)
    }   else if (!is.na(Coculture_df[RowNo, 4]) && !is.na(Coculture_df[RowNo, 10] 
                                                          && abs(Coculture_df[RowNo, 8]) > abs(Coculture_df[RowNo, 14]))) { 
        MatchedPeak_df[RowNo, 1] <- CON2_Name
        MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
        RowNo <- RowNo + 1
        print(RowNo)
    }   else if (is.na(Coculture_df[RowNo, 10]))   {
        MatchedPeak_df[RowNo, 1] <- CON1_Name
        MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 4:9]
        RowNo <- RowNo + 1
        print(RowNo)
    }   else {
        MatchedPeak_df[RowNo, 1] <- CON2_Name
        MatchedPeak_df[RowNo, 2:7] <- Coculture_df[RowNo, 10:15]
        RowNo <- RowNo + 1
        print(RowNo)
    }
}

MatchedPeak_df <- MatchedPeak_df[1:nrow(MatchedPeak_df), ]
df <- setNames(data.frame(matrix(ncol = 1, nrow = nrow(Coculture_df))), "Sample_Ref")
df[1:nrow(df), ] <- Coculture_Name
Refined_Coculture_df <- Coculture_df[1:nrow(Coculture_df), 1:3] #The df containing Coculture data to be cbinded onto.
Refined_Coculture_df <- cbind(df, Refined_Coculture_df)
Refined_Coculture_df <- cbind(Refined_Coculture_df, MatchedPeak_df)
library(dplyr)
library(tidyr)
library(readxl)

#Read_Excel Function to create 3 data frames for the 3 samples to be compared.
Read_Excel <- function(Excel_Name) {
    df_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", Excel_Name, ".xlsm"), 
                          skip = 3)
    
#The following extracts and sorts the top 5 UV maxima from each peak
    UV_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("UV1", "UV2", "UV3", "UV4", "UV5"))
    
    n <- nrow(df_Name)
    i <- 1
    while (i <= n) {
        if (!is.na(df_Name$'UV Peaks')) {
            UV_separate <-df_Name$`UV Peaks`[[i]] ##Looking at row by row. starting with compound 1
            UV_separate <- gsub("\\(","", UV_separate)
            UV_separate <- gsub("\\)","", UV_separate)
            UV_separate <- gsub("< 190","", UV_separate)
            UV_separate <- gsub("s","", UV_separate)
            UV_separate <- strsplit(UV_separate, "\r\n")
            UV_separate <- as.data.frame(UV_separate, col.names = "UV1")
            UV_separate <- separate(UV_separate, 'UV1', c("UV1", "UV1_percent"), sep = " ")
            UV_separate <- UV_separate[order(UV_separate$UV1_percent, decreasing =TRUE),]
            UV_df <- rbind(UV_df, UV_separate[1:5,1])
            i <- i +1    
        }   else {
            i <- i + 1
            UV_separate <- c(NA, NA, NA, NA, NA)
            UV_df <- rbind(UV_df, UV_separate[1:5,1])
        }
    }
    names(UV_df) <- c("UV1", "UV2", "UV3", "UV4", "UV5")   
    UV_df <- transform(UV_df, UV1 = as.numeric(UV1))
    UV_df <- transform(UV_df, UV2 = as.numeric(UV2))
    UV_df <- transform(UV_df, UV3 = as.numeric(UV3))
    UV_df <- transform(UV_df, UV4 = as.numeric(UV4))
    UV_df <- transform(UV_df, UV5 = as.numeric(UV5))
    df_Name <- cbind(df_Name, UV_df)
    df_Name <- select(df_Name, 1:4, 20:24)
    return(df_Name)
}

#Read_UV Function
Read_UV <- function(Excel_Name) {
    UV_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", Excel_Name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(UV_Name)
}
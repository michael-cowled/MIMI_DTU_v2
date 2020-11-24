library(dplyr)
library(tidyr)
library(readxl)

#Read_Excel Function to create 3 data frames for the 3 samples to be compared.
Read_Excel <- function(Excel_Name) {
    df_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", Excel_Name, ".xlsm"), 
                          skip = 3)
    df_Name$`UV Peaks` <- gsub("\\s*\\([^\\)]+\\)","",df_Name$`UV Peaks`)
    df_Name <- separate(df_Name, 'UV Peaks', c("UV1", "UV2", "UV3", "UV4", "UV5"), sep = "\r\n") %>%
        select(1:4, 11:15)
    df_Name <- transform(df_Name, UV1 = as.numeric(UV1))
    df_Name <- transform(df_Name, UV2 = as.numeric(UV2))
    df_Name <- transform(df_Name, UV3 = as.numeric(UV3))
    df_Name <- transform(df_Name, UV4 = as.numeric(UV4))
    df_Name <- transform(df_Name, UV5 = as.numeric(UV5))
    return(df_Name)
}

CON1 <- Read_Excel(CON1_Name)
CON2 <- Read_Excel(CON2_Name)
Coculture <- Read_Excel(Coculture_Name)

#Read_UV Function
Read_UV <- function(Excel_Name) {
    UV_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", Excel_Name, ".xlsm"), sheet = "NormalisedUVVisData")
    return(UV_Name)
}
#NEw testing:
df_Name <- read_excel(paste0("Testing Broad-Scale Interactions/NovaCfiles/", "F1CON", ".xlsm"), 
                      skip = 3)
UV_df <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("UV1", "UV2", "UV3", "UV4", "UV5"))

Matrix_Row_No <- Matrix_Row_No + 1
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
print(UV_df)












library(readxl)

ReadExcel <- function(excel.name) {
    raw.df <- read_excel((excel.name), skip = 3)
    n <- nrow(raw.df)
    new.df <- raw.df[1:n, 1:3]
    new.df2 <- raw.df[1:n, 11]
    raw.data <- cbind(new.df, new.df2)
    return(raw.data)
    
}

ReadUV <- function(excel.name) {
    uv.data <- read_excel((excel.name), sheet = "NormalisedUVVisData")
    return(uv.data)
}

Exporter <- function(FolderName) {
    
    filenames <- list.files(path = paste0("Testing Broad-Scale Interactions/",
                                          FolderName), pattern = "", full.names = TRUE)
    for (i in filenames[1:length(filenames)]) {
        print(i)
        raw.data <- ReadExcel(i)
        uv.data <- ReadUV(i)
        new.name <- gsub(".xlsm","", i)
        write.csv(raw.data, paste0("Simplified/", new.name, "_raw.CSV"), row.names = FALSE)
        write.csv(uv.data, paste0("Simplified/", new.name, "_uv.CSV"), row.names = FALSE)
    }
}

Exporter("NovaCfiles")
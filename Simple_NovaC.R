Exporter <- function(FolderName) {
    
    filenames <- list.files(path = paste0("Testing Broad-Scale Interactions/",
                                          FolderName), pattern = "", full.names = TRUE)
    for (i in filenames[1:length(filenames)]) {
        print(i)
        raw.df <- read_excel(i, skip = 3)
        uv.df <- read_excel(i, sheet = "NormalisedUVVisData")
        new.df <- raw.df[1:nrow(raw.df), 1:3]
        uv.data <- ReadUV(i)
        new.name <- gsub("Testing Broad-Scale Interactions/NovaCfiles/","", i)
        new.name <- gsub(".xlsm","", new.name)
        write.csv(new.df, paste0("Simplified/", new.name, "_raw.CSV"), row.names = FALSE)
        write.csv(uv.data, paste0("Simplified/", new.name, "_uv.CSV"), row.names = FALSE)
    }
}

Exporter("NovaCfiles")
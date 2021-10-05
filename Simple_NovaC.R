## Note: uses the old Read_Excel
library(readxl)

Exporter <- function() {
    
    filenames <- list.files(path = paste0("Testing Broad-Scale Interactions/NovaCfiles/"), 
                            pattern = "", full.names = TRUE)
    print(filenames)
    for (i in filenames[1:length(filenames)]) {
        print(i)
        
        raw.df <- read_excel(i, skip = 3)
        uv.data <- read_excel(i, sheet = "NormalisedUVVisData")
        new.df <- raw.df[1:nrow(raw.df), 1:3]
        
        new.name <- gsub("Testing Broad-Scale Interactions/NovaCfiles/","", i)
        new.name <- gsub(".xlsm","", new.name)
        write.csv(new.df, paste0("Simplified/", new.name, "_raw.CSV"), row.names = FALSE)
        write.csv(uv.data, paste0("Simplified/", new.name, "_uv.CSV"), row.names = FALSE)
    }
    
}

Exporter()
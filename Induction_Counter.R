library(data.table)
library(dplyr)
library(tidyr)

filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/NT_FvA", 
                        pattern = "NT", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)

#For changing NT to F in NT_FvA
full.list$Sample_Ref <- gsub("NT", "F", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("F17", "B1", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("F18", "B2", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("F19", "B3", full.list$Sample_Ref, perl = TRUE)

refined.list <- filter(full.list, Metabolite_Effect == 6)
cocultures <- c(unique(refined.list$Sample_Ref))

induction.df <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in cocultures[1:length(cocultures)]) {
    temp <- filter(refined.list, Sample_Ref == i)
    temp.5perc <- filter(temp, PercArea > 5)
    induction.df <- rbind(induction.df, 
                      c(Sample_Ref = i, Num_Inductions = nrow(temp), 
                        Num_Sig_Inductions = nrow(temp.5perc)))
}

names(induction.df) <- c("Sample_Ref", "Num_Inductions", "Num_Sig_Inductions")
induction.df <- separate(induction.df, Sample_Ref, 
                         into = c("Ref_Culture", "Int_Culture"), sep = "v")
induction.df$Num_Inductions <- as.numeric(induction.df$Num_Inductions)
induction.df$Num_Sig_Inductions <- as.numeric(induction.df$Num_Sig_Inductions)

#Inductions per ref culture

fungi <- c(unique(induction.df$Ref_Culture))
induction.df2 <- data.frame(matrix(ncol = 3, nrow = 0))

for (i in fungi[1:length(fungi)]) {
    temp <- filter(induction.df, Ref_Culture == i)
    induction.df2 <- rbind(induction.df2, c(Ref_Culture = i, 
                            Num_Inductions = sum(temp$Num_Inductions),
                            Num_Sig_Inductions = sum(temp$Num_Sig_Inductions)))
}

names(induction.df2) <- c("Sample_Ref", "Num_Inductions", "Num_Sig_Inductions")
print(induction.df2)
induction.sum <- sum(as.numeric(induction.df2$Num_Inductions))
induction.mean <- induction.sum/length(filenames)
induction.sum
induction.mean

library(data.table)
library(dplyr)
library(tidyr)

filenames <- list.files(path = "Simplified/Testing Broad-Scale Interactions/OutputFiles/NT_Mix", 
                        pattern = "", full.names = TRUE)
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
error <- qnorm(.95)*(sd(induction.df2$Num_Inductions)/sqrt(length(induction.df2$Ref_Culture)))
induction.sum
induction.mean
error
paste0("n is ", length(induction.df2$Ref_Culture))

### Post-use of induction counter a spreadsheet was generated to create a barchart
# Import as 'Inductons' excel object

Inductions <- read_excel("Number_of_Inductions.xlsx")

library(ggplot2)
# Default bar plot
p<- ggplot(Inductions, aes(x=Type, y=Number, fill=Source)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Number-Error, ymax=Number+Error), width=.2,
                  position=position_dodge(.9)) 
print(p)
# Finished bar plot
p+labs(x="Interaction Type", y = "Average Number of Inductions")+
    theme_classic() +
    scale_fill_manual(values=c('#999999','#E69F00'))
##500 x 300
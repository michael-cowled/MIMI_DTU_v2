library(data.table)
library(dplyr)
library(ggplot2)
library(readxl)

mean.list <- setnames(data.frame(matrix(ncol = 8, nrow = 15)), c("mE5", "mE6", "mE5&6", "mAll",
                                                                 "eE5", "eE6", "eE5&6", "eAll"))

for (i in 1:15) {
x <- paste0("NT", i, "v")
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/NT_FvA", pattern = x, full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list <- transform(full.list, Metabolite_Effect = factor(Metabolite_Effect))
full.list <- filter(full.list, Metabolite_Effect != 1)

Effect_5 <- filter(full.list, Metabolite_Effect == 5)
Effect_6 <- filter(full.list, Metabolite_Effect == 6)
Effect_5_6 <- filter(full.list, as.numeric(Metabolite_Effect) > 4)

mean.list[i, 1] <- mean(Effect_5$RetTime_CC)
mean.list[i, 2] <- mean(Effect_6$RetTime_CC)
mean.list[i, 3] <- mean(Effect_5_6$RetTime_CC)
mean.list[i, 4] <- mean(full.list$RetTime_CC)
mean.list[i, 5] <- error <- qnorm(.95)*(sd(Effect_5$RetTime_CC)/sqrt(length(Effect_5$RetTime_CC)))
mean.list[i, 6] <- error <- qnorm(.95)*(sd(Effect_6$RetTime_CC)/sqrt(length(Effect_6$RetTime_CC)))
mean.list[i, 7] <- error <- qnorm(.95)*(sd(Effect_5_6$RetTime_CC)/sqrt(length(Effect_5_6$RetTime_CC)))
mean.list[i, 8] <- error <- qnorm(.95)*(sd(full.list$RetTime_CC)/sqrt(length(full.list$RetTime_CC)))

}

for (i in 1:15) {
    col1 <- c(paste0("F", i), paste0("F", i), paste0("F", i), paste0("F", i))
    if (i > 1) {
        extraction2 <- cbind(col1, c("E5", "E6", "E5 & E6", "All"), 
                                              transpose(mean.list[i,1:4]), transpose(mean.list[i,5:8]))
        names(extraction2) <- names(extraction)
        extraction <- rbind(extraction, extraction2)
    } else {
        extraction <- cbind(col1, c("E5", "E6", "E5 & E6", "All"), 
                            transpose(mean.list[i,1:4]), transpose(mean.list[i,5:8]))
        names(extraction) <- c("fungus", "Effect", "mean", "se")
        }
}

extraction$fungus <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", extraction$fungus, perl = TRUE)
Fungi_Assignments <- read_excel("Testing Broad-Scale Interactions/Fungi_Assignments_nt4.xlsx")
extraction <- merge(extraction, Fungi_Assignments, by="fungus", all.x=TRUE)

p<- ggplot(extraction, aes(x=Ref_Species, y=mean, fill=Effect)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4,
                  position=position_dodge(.9)) 
barchart <- p+labs(x="Reference Fungus", y = "Retention Time")+
    theme_classic() +
    theme(axis.text.x=element_text(angle=45, hjust = 1))
barchart
ggsave(barchart,filename="barchart_fvf_tal.png",height=3.5,width=10.20,units="in",dpi=200)
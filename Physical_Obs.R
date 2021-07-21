# Based on part 1 of interaction report
# Import Physical_Obs object

ratio.df3[order(ratio.df3$Int_Culture),]

Physical_Obs$Ref_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", Physical_Obs$Ref_Culture, perl = TRUE)
Physical_Obs$Int_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", Physical_Obs$Int_Culture, perl = TRUE)
ratio.df4 <- mutate(ratio.df3, Physical_Effect = Physical_Obs$Physical_Effect)
    
inhibited <- filter(ratio.df4, Physical_Effect == 1)
inhibiting <- filter(ratio.df4, Physical_Effect == 2)
contact.inhibition <- filter(ratio.df4, Physical_Effect == "C")
phenotype <- filter(ratio.df4, Physical_Effect == "P")
no.change <- filter(ratio.df4, Physical_Effect == "N")

mean(inhibited$Avg_Ratio)
mean(inhibiting$Avg_Ratio)
mean(contact.inhibition$Avg_Ratio)
mean(phenotype$Avg_Ratio)
mean(no.change$Avg_Ratio)
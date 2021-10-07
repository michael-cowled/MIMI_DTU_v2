# Based on part 1 of interaction report
# Import Physical_Obs object

ratio.df3[order(ratio.df3$Int_Culture),]

Physical_Obs$Ref_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", Physical_Obs$Ref_Culture, perl = TRUE)
Physical_Obs$Int_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", Physical_Obs$Int_Culture, perl = TRUE)
ratio.df4 <- mutate(ratio.df3, Physical_Effect = Physical_Obs$Physical_Effect)
    
inhibited <- filter(ratio.df4, Physical_Effect == 1) %>%
    mutate(Norm_Avg_Ratio = Avg_Ratio - mean(ratio.df4$Avg_Ratio))
inhibiting <- filter(ratio.df4, Physical_Effect == 2)  %>%
    mutate(Norm_Avg_Ratio = Avg_Ratio - mean(ratio.df4$Avg_Ratio))
contact.inhibition <- filter(ratio.df4, Physical_Effect == "C") %>%
    mutate(Norm_Avg_Ratio = Avg_Ratio - mean(ratio.df4$Avg_Ratio))
phenotype <- filter(ratio.df4, Physical_Effect == "P") %>%
    mutate(Norm_Avg_Ratio = Avg_Ratio - mean(ratio.df4$Avg_Ratio))
no.change <- filter(ratio.df4, Physical_Effect == "N") %>%
    mutate(Norm_Avg_Ratio = Avg_Ratio - mean(ratio.df4$Avg_Ratio))

mean(inhibited$Norm_Avg_Ratio)
mean(inhibiting$Norm_Avg_Ratio)
mean(contact.inhibition$Norm_Avg_Ratio)
mean(phenotype$Norm_Avg_Ratio)
mean(no.change$Norm_Avg_Ratio)

res <- t.test(phenotype$Norm_Avg_Ratio, no.change$Norm_Avg_Ratio, alternative = "greater", var.equal = FALSE)
print(res$p.value)

qnorm(.95)*(sd(inhibited$Norm_Avg_Ratio)/sqrt(length(inhibited$Norm_Avg_Ratio)))
qnorm(.95)*(sd(inhibiting$Norm_Avg_Ratio)/sqrt(length(inhibiting$Norm_Avg_Ratio)))
qnorm(.95)*(sd(contact.inhibition$Norm_Avg_Ratio)/sqrt(length(contact.inhibition$Norm_Avg_Ratio)))
qnorm(.95)*(sd(phenotype$Norm_Avg_Ratio)/sqrt(length(phenotype$Norm_Avg_Ratio)))
qnorm(.95)*(sd(no.change$Norm_Avg_Ratio)/sqrt(length(no.change$Norm_Avg_Ratio)))

## Creation of barchart
## Import an object with the data

library(ggplot2)
# Default bar plot
p<- ggplot(Mean_Avg_Ratio_NT, aes(x=Phys_Obs, y=Number)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=Number-Error, ymax=Number+Error), width=.2,
                  position=position_dodge(.9)) 
print(p)
# Finished bar plot
p+labs(x="Interaction Type", y = "Average Number of Inductions")+
    theme_classic() +
    scale_fill_manual(values=c('#999999','#E69F00'))
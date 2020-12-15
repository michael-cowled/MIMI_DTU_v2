library(data.table)
library(ggplot2)

#Combines all output files into a list of all metabolites from all cocultures

filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/", pattern = "F", full.names = TRUE)
my_data <- lapply(filenames, read.csv)
full_list <- rbindlist(my_data, use.names=TRUE, fill=FALSE)
full_list <- transform(full_list, Metabolite_Effect = factor(Metabolite_Effect))

#Filters into different lists based on metabolite_effect
Effect_1 <- filter(full_list, Metabolite_Effect == 1)
Effect_2 <- filter(full_list, Metabolite_Effect == 2)
Effect_3 <- filter(full_list, Metabolite_Effect == 3)
Effect_4 <- filter(full_list, Metabolite_Effect == 4)
Effect_5 <- filter(full_list, Metabolite_Effect == 5)
Effect_6 <- filter(full_list, Metabolite_Effect == 6)

#Generates a sequence of histograms comparing metabolite effect to ret_time
par(mfrow=c(2,3))
hist(Effect_1$RetTime_CON, main = "Complete Suppression")
hist(Effect_2$RetTime_CC, main = "Suppression")
hist(Effect_3$RetTime_CC, main = "Little to No Change")
hist(Effect_4$RetTime_CC, main = "Enhancement")
hist(Effect_5$RetTime_CC, main = "Major Enhancement")
hist(Effect_6$RetTime_CC, main = "Induction or Unmatched")

par(mfrow=c(1,1))
a <- hist(Effect_6$RetTime_CC, main = "Induction or Unmatched")
b <- hist(full_list$RetTime_CC, main = "Induction or Unmatched")
plot( a, col=rgb(0,0,1,1/4), xlim=c(0,11), ylim=c(0,800), main = "Overlay of Induction to All Metabolites")  # first histogram
plot( b, col=rgb(1,0,0,1/4), xlim=c(0,11), ylim=c(0,800), add=T)  # second
plot( c, col=rgb(0,1,0,1/4), xlim=c(0,11), ylim=c(0,800), add=T)  # third


#Separate lists based on Coculturing Fungus
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/", pattern = "F1v", full.names = TRUE)
my_data <- lapply(filenames, read.csv)
F1v <- rbindlist(my_data, use.names=TRUE, fill=FALSE)
F1v <- transform(F1v, Metabolite_Effect = factor(Metabolite_Effect))

#e.g. boxplot
boxplot(full_list$RetTime_CON ~ full_list$Metabolite_Effect)


#Separate lists based on Coculturing Fungus
par(mfrow=c(3,5), mar=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
for (i in 1:15) {
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/", pattern = paste0("F", i, "v"), full.names = TRUE)
my_data <- lapply(filenames, read.csv)
a <- rbindlist(my_data, use.names=TRUE, fill=FALSE)
a  <- transform(a , Metabolite_Effect = factor(Metabolite_Effect))
boxplot(a$RetTime_CC ~ a$Metabolite_Effect, main = paste0("F", i, "v"), xlab="", ylab="") }


#SCatterplot with colour
full_list <- filter(full_list, PeakRatio != -100)
full_list <- transform(full_list, Sample_Ref = factor(Sample_Ref))
full_list <- transform(full_list, Matched_CON = factor(Matched_CON))
qplot(full_list$RetTime_CC, full_list$PeakRatio, color = full_list$Sample_Ref)
qplot(full_list$RetTime_CC, full_list$PeakRatio, color = full_list$Matched_CON)

a <- filter(full_list, PeakRatio >2500)
a
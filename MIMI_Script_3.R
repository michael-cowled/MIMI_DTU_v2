##########################################################
##########################################################
# Microbial Interaction Metabolite Integrator - Script 3 #
##########################################################
##########################################################

# COPYWRIGHT: © Macquarie University - Michael Cowled and contributors [2021]

# INPUT: The output files from Script_2.

# OUTPUT: Graphs and plots either in report fashion (as indicated by Script_3),
# or generated from individual functions based on the needs of the user.


# R packages required to be loaded in:

library(data.table)
library(ggplot2)

#------------------------------------------------------------------------------#

#############################################
# 1.ReadOutput
#############################################

# Combines all output files into a list of all metabolites from all cocultures

# Args:
    # dataset - specifies a subsetted dataset folder. Leave blank if N/A.

# Returns:
    # full.list - a dataframe consisting of ALL output files combined into one.

ReadOutput <- function(dataset) {
filenames <- list.files(path = paste0("Testing Broad-Scale Interactions/OutputFiles/", dataset),
                        pattern = "F", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list <- transform(full.list, Metabolite_Effect = factor(Metabolite_Effect))
return(full.list)
}

# Subsets full.list into specific effect categories
effect.2 <- all.other.effects <- filter(full.list, as.numeric(Metabolite_Effect) < 3 
                                        & as.numeric(Metabolite_Effect) > 1)
effect.5.and.6 <- filter(full.list, as.numeric(Metabolite_Effect) > 4)
all.other.effects <- filter(full.list, as.numeric(Metabolite_Effect) < 5 
                            & as.numeric(Metabolite_Effect) > 2)

# Generates a side-by-side histogram comparison of RetTime vs. # of metabolites

par(mfrow=c(1,3))
hist(effect.2$RetTime_CC, main = "Suppressed", 
     xlab = "Retention Time", ylab = "Number of Secondary Metabolites")
hist(effect.5.and.6$RetTime_CC, main = "Induction or Major Enhancement", 
     xlab = "Retention Time", ylab = "Number of Secondary Metabolites")
hist(all.other.effects$RetTime_CC, main = "All other effects", 
     xlab = "Retention Time", ylab = "Number of Secondary Metabolites")

# Produces a df comparing the number of metabolites in specific RetTime ranges.

a <- nrow(filter(effect.5.and.6, RetTime_CC >= 2 & RetTime_CC <= 7))
b <- nrow(effect.5.and.6) - a
c <- nrow(filter(all.other.effects, RetTime_CC >= 2 & RetTime_CC <= 7))
d <- nrow(all.other.effects) - c
Metabolite_Count <- c(a, b, c, d)
Range <- c("2-7 min", "0-2 min, 7-11 min", "2-7 min", "0-2 min, 7-11 min")
Condition <- c("Induced or majorly enhanced", "Induced or majorly enhanced", 
               "Other effects", "Other effects")
hist.data <- data.frame(Metabolite_Count, Range, Condition)

# Generates the corresponding histogram

ggplot(hist.data, aes(fill=Range, y=Metabolite_Count, x=Condition)) + 
    geom_bar(position="fill", stat="identity") + xlab("Effect") + 
    ylab("Proportion of Secondary Metabolites") + 
    scale_fill_discrete(name = "Retention Time Range") + theme_minimal() + 
    theme(legend.position = "bottom")

# Comparing Inductions to all other effects in 3-4 min range

effect.6 <- filter(full.list, as.numeric(Metabolite_Effect) > 5)
all.other.effects <- filter(full.list, as.numeric(Metabolite_Effect) < 6 
                            & as.numeric(Metabolite_Effect) > 1)
a <- nrow(filter(effect.6, RetTime_CC >= 3 & RetTime_CC <= 4))
b <- nrow(effect.6)
c <- nrow(filter(all.other.effects, RetTime_CC >= 3 & RetTime_CC <= 4))
d <- nrow(all.other.effects)
Metabolite_Count <- c(a, b, c, d)
Range <- c("3-4 min", "0-3 min, 4-11 min", "3-4 min", "0-3 min, 4-11 min")
Condition <- c("Induced", "Induced", 
               "Other effects", "Other effects")
hist.data <- data.frame(Metabolite_Count, Range, Condition)

# Boxplot comparisons of metabolite effect vs retention time

full.list.minus.effect.1 <- filter(full.list, as.numeric(Metabolite_Effect) <7 & as.numeric(Metabolite_Effect) > 1)
effect.names <- c("Suppression", "Little to no change", "Enhancement", "Major enhancement", "Induction")
ggplot(data=full.list.minus.effect.1, aes(x=Metabolite_Effect, y=RetTime_CC)) + geom_boxplot(aes(fill=Metabolite_Effect)) + 
    ylab("Retention Time (min)") + xlab("Metabolite Effect") +
    stat_summary(fun.y=mean, geom="point", shape=1, size=3) +
    theme(axis.text.x = element_text(angle = 90, face = "italic")) + scale_x_discrete(labels = effect.names)
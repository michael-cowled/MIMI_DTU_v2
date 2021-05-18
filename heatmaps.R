library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/Tal_AvA", pattern = "", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)

#For changing NT to F in NT_FvA
full.list$Sample_Ref <- gsub("NT", "F", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("F17", "B1", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("F18", "B2", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("F19", "B3", full.list$Sample_Ref, perl = TRUE)
#Correcting for culture "86280v3"
full.list[full.list == "86280v3vRA14283c"] <- "86280V3vRA14283c"
full.list[full.list == "86280v3vMA9095"] <- "86280V3vMA9095"
full.list[full.list == "86280v3vAS6166"] <- "86280V3vAS6166"
full.list[full.list == "86280v3vAS5549"] <- "86280V3vAS5549"
full.list[full.list == "86280v3vAS4461"] <- "86280V3vAS4461"
full.list[full.list == "86280v3v8651"] <- "86280V3v8651"

refined.list <- filter(full.list, Metabolite_Effect < 6 & Metabolite_Effect > 1)
Cocultures <- c(unique(refined.list$Sample_Ref))

ratio.df <- data.frame(matrix(ncol = 3, nrow = 0))
for (i in Cocultures[1:length(Cocultures)]) {
    print(i)
    temp <- filter(refined.list, Sample_Ref == i)
    ratio.df <- rbind(ratio.df, 
                      c(Sample_Ref = i, Max_Ratio = max(temp$PeakRatio), 
                        Min_Ratio = min(temp$PeakRatio)))
}

names(ratio.df) <- c("Sample_Ref", "Max_Ratio", "Min_Ratio")
ratio.df <- separate(ratio.df, Sample_Ref, into = c("Ref_Culture", "Int_Culture"), sep = "v")
ratio.df$Max_Ratio <- as.numeric(ratio.df$Max_Ratio)
ratio.df$Min_Ratio <- as.numeric(ratio.df$Min_Ratio)

ratio.df <- mutate(ratio.df, Log_Max_Ratio = log(Max_Ratio))

#For NT fungi
ratio.df$Ref_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", ratio.df$Ref_Culture, perl = TRUE)
ratio.df$Int_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", ratio.df$Int_Culture, perl = TRUE)

heatmap.fvf.enhancements <- ggplot(data = ratio.df, aes(x=Int_Culture, y=Ref_Culture)) + 
    geom_tile(aes(fill=Max_Ratio)) +
    scale_fill_gradient2(low = "purple", high = "dark red", mid = "white", 
                         midpoint = min(ratio.df$Max_Ratio), space = "Lab", name="Ratio (%)") +
    labs(x="Interacting Culture", y="Reference Culture") +
    theme_grey(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust = 1))
heatmap.fvf.enhancements
ggsave(heatmap.fvf.enhancements,filename="heatmap.fvf.enhancements.png",height=2.24,width=3,units="in",dpi=200)

heatmap.fvf.log.enhancements <- ggplot(data = ratio.df, aes(x=Int_Culture, y=Ref_Culture)) + 
    geom_tile(aes(fill=Log_Max_Ratio)) +
    scale_fill_gradient2(low = "white", high = "dark red", mid = "white", 
                         midpoint = 0, space = "Lab", name="log Ratio", limits = c(0, 10)) +
    labs(x="Interacting Culture", y="Reference Culture") +
    theme_grey(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust = 1))
heatmap.fvf.enhancements
ggsave(heatmap.fvf.log.enhancements,filename="heatmap.fvf.log.enhancements.png",height=2.24,width=3,units="in",dpi=200)

heatmap.fvf.suppressions <- ggplot(data = ratio.df, aes(x=Int_Culture, y=Ref_Culture)) + 
    geom_tile(aes(fill=Min_Ratio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = max(ratio.df$Min_Ratio), space = "Lab", name="Ratio (%)") +
    labs(x="Interacting Culture", y="Reference Culture") +
    theme_grey(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust = 1))
heatmap.fvf.suppressions
ggsave(heatmap.fvf.suppressions,filename="heatmap.fvf.suppressions.png",height=2.24,width=3,units="in",dpi=200)

#Avg Ratio Heatmap

refined.list <- filter(full.list, Metabolite_Effect < 5 & Metabolite_Effect > 1)
cond.5.list <- filter(full.list, Metabolite_Effect == 5)
cond.5.list$PeakRatio <- 100
refined.list <- rbind(refined.list, cond.5.list)
Cocultures <- c(unique(refined.list$Sample_Ref))

ratio.df <- data.frame(matrix(ncol = 2, nrow = 0))
for (i in Cocultures[1:length(Cocultures)]) {
    print(i)
    temp <- filter(refined.list, Sample_Ref == i)
    ratio.df <- rbind(ratio.df, 
                      c(Sample_Ref = i, Avg_Ratio = mean(temp$PeakRatio)))
}

names(ratio.df) <- c("Sample_Ref", "Avg_Ratio")
ratio.df <- separate(ratio.df, Sample_Ref, into = c("Ref_Culture", "Int_Culture"), sep = "v")
ratio.df$Avg_Ratio <- as.numeric(ratio.df$Avg_Ratio)

#For NT fungi
ratio.df$Ref_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", ratio.df$Ref_Culture, perl = TRUE)
ratio.df$Int_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", ratio.df$Int_Culture, perl = TRUE)

#Correcting for "86280v3"
ratio.df[ratio.df == "86280V03"] <- "82680"

heatmap.fvf.avg <- ggplot(data = ratio.df, aes(x=Int_Culture, y=Ref_Culture)) + 
    geom_tile(aes(fill=Avg_Ratio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = 0, space = "Lab", name="Ratio (%)") +
    labs(x="Interacting Culture", y="Reference Culture") +
    theme_grey(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust = 1))
heatmap.fvf.avg
ggsave(heatmap.fvf.avg,filename="heatmap.fvf.avg.png",height=2.24,width=3,units="in",dpi=200)

mean(ratio.df$Avg_Ratio, na.rm = TRUE)
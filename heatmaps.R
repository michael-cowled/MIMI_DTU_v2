filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/Tal_FvF", pattern = "F", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list <- transform(full.list, Metabolite_Effect = factor(Metabolite_Effect))
ALEI <- filter(full.list, as.numeric(Metabolite_Effect) < 5 & as.numeric(Metabolite_Effect) > 1)
Cocultures <- c(unique(ALEI$Sample_Ref))

ratio.df <- data.frame(matrix(ncol = 4, nrow = 0))
for (i in Cocultures[1:length(Cocultures)]) {
    print(i)
    temp <- filter(ALEI, Sample_Ref == i)
    ratio.df <- rbind(ratio.df, 
                      c(Sample_Ref = i, Max_Ratio = max(temp$PeakRatio), 
                        Min_Ratio = min(temp$PeakRatio), Avg_Ratio = mean(temp$PeakRatio)))
}

names(ratio.df) <- c("Sample_Ref", "Max_Ratio", "Min_Ratio", "Avg_Ratio")
ratio.df <- separate(ratio.df, Sample_Ref, into = c("Ref_Culture", "Int_Culture"), sep = "v")
ratio.df$Max_Ratio <- as.numeric(ratio.df$Max_Ratio)
ratio.df$Min_Ratio <- as.numeric(ratio.df$Min_Ratio)
ratio.df$Avg_Ratio <- as.numeric(ratio.df$Avg_Ratio)

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

heatmap.fvf.suppressions <- ggplot(data = ratio.df, aes(x=Int_Culture, y=Ref_Culture)) + 
    geom_tile(aes(fill=Min_Ratio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = max(ratio.df$Min_Ratio), space = "Lab", name="Ratio (%)") +
    labs(x="Interacting Culture", y="Reference Culture") +
    theme_grey(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust = 1))
heatmap.fvf.suppressions
ggsave(heatmap.fvf.suppressions,filename="heatmap.fvf.suppressions.png",height=2.24,width=3,units="in",dpi=200)

heatmap.fvf.avg <- ggplot(data = ratio.df, aes(x=Int_Culture, y=Ref_Culture)) + 
    geom_tile(aes(fill=Avg_Ratio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = 0, space = "Lab", name="Ratio (%)") +
    labs(x="Interacting Culture", y="Reference Culture") +
    theme_grey(base_size=8) +
    theme(axis.text.x=element_text(angle=45, hjust = 1))
heatmap.fvf.avg
ggsave(heatmap.fvf.avg,filename="heatmap.fvf.avg.png",height=2.24,width=3,units="in",dpi=200)
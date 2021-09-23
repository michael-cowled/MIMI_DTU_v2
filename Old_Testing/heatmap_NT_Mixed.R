filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/NT_fvf_mixed", 
                        pattern = "F1", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list$Sample_Ref <- gsub("F1v", "", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", full.list$Sample_Ref, perl = TRUE)
subset.list <- filter(full.list, Matched_con == "F1CON") %>%
    select(Sample_Ref, PeakNo_con, RetTime_con, Matched_con, PeakRatio) %>%
    filter(Sample_Ref != "F03" & Sample_Ref != "F05" & Sample_Ref != "F14") %>%
    mutate(logPeakRatio = log((PeakRatio/100) + 1.01))

heatmap_fvf_mixed <- ggplot(data = subset.list, aes(x=PeakNo_con, y=Sample_Ref)) + 
    geom_tile(aes(fill=logPeakRatio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = 0, limit = c(min(subset.list$logPeakRatio), 
                                                 max(subset.list$logPeakRatio)), 
                         space = "Lab", name="log(%PeakArea)") +
    labs(x="Peak Number", y="Interacting Culture") +
    theme_grey(base_size=8)
heatmap_fvf_mixed
ggsave(heatmap_fvf_mixed,filename="heatmap_fvf_mixed.png",height=2.24,width=5,units="in",dpi=200)

# Making an accompany table to go alongside heatmap

unique.peaks <- unique(subset.list$PeakNo_con)
unique.peaks <- sort(unique.peaks)
summary.table <- data.frame(matrix(ncol = 2, nrow = length(unique.peaks)))
summary.table[,1] <- unique.peaks

for (i in unique.peaks[1:length(unique.peaks)]) {
    temp <- filter(subset.list, PeakNo_con == i)
    tempRT <- temp$RetTime_con[1]
    summary.table[i,2] <- tempRT
}
summary.table <- setNames(summary.table, c("PeakNo", "RetTime"))
print(summary.table, digits = 2)


###testing (clustered heatmap)
df.as.matrix <- select(subset.list, Sample_Ref, PeakNo_con, logPeakRatio) 
df.as.matrix <- spread(df.as.matrix, PeakNo_con, logPeakRatio)
Fungi_Assignments <- read_excel("Testing Broad-Scale Interactions/Fungi_Assignments_nt4.xlsx")
df.as.matrix <- merge(df.as.matrix, Fungi_Assignments, by="Sample_Ref", all.x=TRUE)
names <- as.vector(df.as.matrix$Species)
df.as.matrix <- select(df.as.matrix, -Sample_Ref, -Species)

df.as.matrix <- as.matrix(df.as.matrix)
row.names(df.as.matrix) <- names

colfunc <- colorRampPalette(c("dark blue", "white", "dark red"))

heatmap.2(df.as.matrix, col=colfunc(15), xlab = "Peak Number", 
          scale = "none", trace = "none", key=TRUE , key.ylab="Count", 
          key.xlab = expression(paste("log(", alpha, ")")), key.title="Colour Key & Histogram",
          key.ytickfun = FALSE, density.info = "none")
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/NT_fvf_mixed", 
                        pattern = "F1", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list$Sample_Ref <- gsub("F1v", "", full.list$Sample_Ref, perl = TRUE)
full.list$Sample_Ref <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", full.list$Sample_Ref, perl = TRUE)
subset.list <- filter(full.list, Matched_con == "F1CON") %>%
    select(Sample_Ref, PeakNo_con, PeakRatio) %>%
    mutate(logPeakRatio = log((PeakRatio/100) + 1))

heatmap_fvf <- ggplot(data = subset.list, aes(x=PeakNo_con, y=Sample_Ref)) + 
    geom_tile(aes(fill=logPeakRatio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = 0, limit = c(min(subset.list$logPeakRatio), 
                                                 max(subset.list$logPeakRatio)), 
                         space = "Lab", name="log(%PeakArea)") +
    labs(x="Peak Number", y="Interacting Culture") +
    theme_grey(base_size=8)
heatmap_fvf
ggsave(heatmap_fvf,filename="heatmap_fvf.png",height=1.75,width=10.20,units="in",dpi=200)
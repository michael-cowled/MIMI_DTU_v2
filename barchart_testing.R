filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/Tal_FvF", pattern = "F", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list.tal <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list.tal <- transform(full.list.tal, Metabolite_Effect = factor(Metabolite_Effect))

filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/NT_FvF", pattern = "F", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list.nt <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list.nt <- transform(full.list.nt, Metabolite_Effect = factor(Metabolite_Effect))

boxplot(full.list.tal$RetTime_CC ~ full.list$Metabolite_Effect)
boxplot(full.list.tal$RetTime_CC ~ full.list$Metabolite_Effect, las = 2,
        ylab = "Retention Time (min)", xlim = c(1.5,6.5), 
        names = c("Not present", "Suppression", "Little to no change", "Enhancement", "Major enhancement", "Induction"))
mtext("Metabolite Effect", side = 1, line = 7)

# Produce a df comparing the proportion of metabolites in each category
effect <- c(1:6)
tal <- c(nrow(filter(full.list.tal, Metabolite_Effect == 1)), nrow(filter(full.list.tal, Metabolite_Effect == 2)),
         nrow(filter(full.list.tal, Metabolite_Effect == 3)), nrow(filter(full.list.tal, Metabolite_Effect == 4)), 
         nrow(filter(full.list.tal, Metabolite_Effect == 5)), nrow(filter(full.list.tal, Metabolite_Effect == 6)))
nt <- c(nrow(filter(full.list.nt, Metabolite_Effect == 1)), nrow(filter(full.list.nt, Metabolite_Effect == 2)),
        nrow(filter(full.list.nt, Metabolite_Effect == 3)), nrow(filter(full.list.nt, Metabolite_Effect == 4)), 
        nrow(filter(full.list.nt, Metabolite_Effect == 5)), nrow(filter(full.list.nt, Metabolite_Effect == 6)))
df <- data.frame(matrix(ncol = 3, nrow = 6))
df[1:6,1] <- effect
df[1:6,2] <- tal
df[1:6,3] <- nt
names(df) <- c("effect", "talented", "nt")
df <- mutate(df, Tal_prop = talented*100/sum(df$talented), nt_prop = nt*100/sum(df$nt))

#Further dataframe manipulations
group <- c(rep("Talented", 5), rep("NT", 5))
effect2 <- c(rep(effect, 2))
percent <- c(df$Tal_prop, df$nt_prop)
data <- data.frame(group, effect2, percent)

ggplot(data, aes(fill=effect2, y=percent, x=group)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = F) +
    ggtitle("Difference in Effects betweeen Talented and NT Fungal Interactions") +
    theme(axis.text.x = element_text(angle = 0, face = "italic")) + 
    scale_x_discrete(labels = c("Talented", "NT"))
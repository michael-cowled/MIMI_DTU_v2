library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggdendro)
library(reshape2)
library(viridis)
library(RColorBrewer)
library(gplots)

#Combines all output files into a list of all metabolites from all cocultures

filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/NT_FvF", pattern = "F", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list <- transform(full.list, Metabolite_Effect = factor(Metabolite_Effect))

#Filters into different lists based on metabolite_effect
Effect_1 <- filter(full.list, Metabolite_Effect == 1)
Effect_2 <- filter(full.list, Metabolite_Effect == 2)
Effect_3 <- filter(full.list, Metabolite_Effect == 3)
Effect_4 <- filter(full.list, Metabolite_Effect == 4)
Effect_5 <- filter(full.list, Metabolite_Effect == 5)
Effect_6 <- filter(full.list, Metabolite_Effect == 6)
All_Effects_Less_Induction <- filter(full.list, as.numeric(Metabolite_Effect) < 5 & as.numeric(Metabolite_Effect) > 1)
Effect_5_6 <- filter(full.list, as.numeric(Metabolite_Effect) > 4)

#Generates a sequence of histograms comparing metabolite effect to ret_time
par(mfrow=c(2,3))
hist(Effect_1$RetTime_CON, main = "Complete Suppression")
hist(Effect_2$RetTime_CC, main = "Suppression")
hist(Effect_3$RetTime_CC, main = "Little to No Change")
hist(Effect_4$RetTime_CC, main = "Enhancement")
hist(Effect_5$RetTime_CC, main = "Major Enhancement")
hist(Effect_6$RetTime_CC, main = "Induction or Unmatched")
hist(All_Effects_Less_Induction$RetTime_CC, main = "All other effects", xlab = "Retention Time", ylab = "Number of Secondary Metabolites")
hist(Effect_5_6$RetTime_CC, main = "Induction or Major Enhancement", xlab = "Retention Time", ylab = "Number of Secondary Metabolites")

par(mfrow=c(1,1))
a <- hist(Effect_6$RetTime_CC, main = "Induction or Unmatched")
b <- hist(full.list$RetTime_CC, main = "Full List")
c <- hist(All_Effects_Less_Induction$RetTime_CC, main = "All other effects")
d <- hist(Effect_5_6$RetTime_CC, main = "Induction or Unmatched")
plot( a, col=rgb(0,0,1,1/4), xlim=c(0,11), ylim=c(0,600), main = "Overlay of Induction to All Metabolites")  # first histogram
plot( b, col=rgb(1,0,0,1/4), xlim=c(0,11), ylim=c(0,800), add=T)  # second
plot( c, col=rgb(0,1,0,1/4), xlim=c(0,11), ylim=c(0,600), add=T)  # third
plot( d, col=rgb(1,0,1,1/4), xlim=c(0,11), ylim=c(0,600), add=T)  # third

#Separate lists based on Coculturing Fungus
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/", pattern = "F1v", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
F1v <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
F1v <- transform(F1v, Metabolite_Effect = factor(Metabolite_Effect))

#e.g. boxplot
boxplot(full.list$RetTime_CC ~ full.list$Metabolite_Effect)
boxplot(full.list$RetTime_CC ~ full.list$Metabolite_Effect, las = 2,
        ylab = "Retention Time (min)", xlim = c(1.5,6.5), 
        names = c("Not present", "Suppression", "Little to no change", "Enhancement", "Major enhancement", "Induction"))
mtext("Metabolite Effect", side = 1, line = 7)

full.list.minus.effect.1 <- filter(full.list, as.numeric(Metabolite_Effect) <7 & as.numeric(Metabolite_Effect) > 1)
effect.names <- c("Suppression", "Little to no change", "Enhancement", "Major enhancement", "Induction")
ggplot(data=full.list.minus.effect.1, aes(x=Metabolite_Effect, y=RetTime_CC)) + geom_boxplot(aes(fill=Metabolite_Effect)) + 
    ylab("Retention Time (min)") + xlab("Metabolite Effect") +
    stat_summary(fun.y=mean, geom="point", shape=1, size=4) +
    theme(axis.text.x = element_text(angle = 90, face = "italic")) + scale_x_discrete(labels = effect.names)

#Separate lists based on Coculturing Fungus
par(mfrow=c(3,5), mar=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
for (i in 1:15) {
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/", pattern = paste0("F", i, "v"), full.names = TRUE)
my.data <- lapply(filenames, read.csv)
a <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
a  <- transform(a , Metabolite_Effect = factor(Metabolite_Effect))
boxplot(a$RetTime_CC ~ a$Metabolite_Effect, main = paste0("F", i, "v"), xlab="", ylab="") }

#SCatterplot with colour
full.list <- filter(full.list, PeakRatio != -100)
full.list <- transform(full.list, Sample_Ref = factor(Sample_Ref))
full.list <- transform(full.list, Matched_con = factor(Matched_con))
qplot(full.list$RetTime_CC, full.list$PeakRatio, color = full.list$Sample_Ref)
qplot(full.list$RetTime_CC, full.list$PeakRatio, color = full.list$Matched_con)

a <- filter(full.list, PeakRatio >2500)
a


#For looking at a single fungus
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/NT_fvf_mixed", 
                        pattern = "F1", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)

subset.list <- filter(full.list, Matched_con == "F1CON") %>%
    select(Sample_Ref, PeakNo_con, PeakRatio) %>%
    mutate(logPeakRatio = log((PeakRatio/100) + 1))

heatmap_fvf <- ggplot(data = subset.list, aes(x=PeakNo_con, y=Sample_Ref)) + 
    geom_tile(aes(fill=logPeakRatio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = 0, limit = c(min(subset.list$logPeakRatio), 
                                                 max(subset.list$logPeakRatio)), 
                         space = "Lab", name="log(%PeakArea)") +
    labs(x="Peak Number", y="Coculture") +
    theme_grey(base_size=8)
heatmap_fvf
ggsave(heatmap_fvf,filename="heatmap_fvf.png",height=1.75,width=10.20,units="in",dpi=200)

#By countfactor
subset.list <- filter(full.list, Matched_con == "FP1927CON", Metabolite_Effect == 1 | Metabolite_Effect == 6) %>%
    select(Sample_Ref, PeakNo_con, PeakRatio, Metabolite_Effect) %>%
    mutate(logPeakRatio = log((PeakRatio/100) + 1)) %>%
mutate(countfactor=cut(PeakRatio,breaks=c(-100,-10,10,100,max(PeakRatio,na.rm=T)),
                       labels=c("suppression","no change","minor enhancement","major enhancement"))) %>%
    mutate(countfactor=factor(as.character(countfactor),levels=rev(levels(countfactor))))

heatmap_fvf <- ggplot(data = subset.list, aes(x=PeakNo_con, y=Sample_Ref)) + 
    geom_tile(aes(fill=countfactor)) +
    scale_fill_manual(values=c("dark red","dark orange", "white","dark blue"),na.value = "grey90")+
    labs(x="Peak Number", y="Coculture") +
    theme_grey(base_size=8)
heatmap_fvf
ggsave(heatmap_fvf,filename="heatmap_fvf.png",height=5.5,width=8.8,units="in",dpi=200)

#Just inductions
filenames <- list.files(path = "Testing Broad-Scale Interactions/OutputFiles/Tal_fvf_mixed", pattern = "F", full.names = TRUE)
my.data <- lapply(filenames, read.csv)
full.list <- rbindlist(my.data, use.names=TRUE, fill=FALSE)
full.list$Metabolite_Effect <- as.numeric(full.list$Metabolite_Effect)
subset.list <- filter(full.list, is.na(Matched_con) & Metabolite_Effect == 6) %>%
    select(Sample_Ref, PeakNo_CC, PeakRatio, Metabolite_Effect)


heatmap_fvf_inductions <- ggplot(data = subset.list, aes(x=PeakNo_CC, y=Sample_Ref)) + 
    geom_tile(aes(fill=Metabolite_Effect)) +
    scale_fill_gradient2(low = "purple", high = "dark green", mid = "dark green", 
                         midpoint = 0, space = "Lab", name="Inductions") +
    labs(x="Peak Number", y="Coculture") +
    theme_grey(base_size=8)
heatmap_fvf_inductions
ggsave(heatmap_fvf_inductions,filename="heatmap_fvf_inductions.png",height=1.75,width=10.20,units="in",dpi=200)
    

## Testing for single fungus loop
# Post importing full.list

coculture.list <- filter(full.list, is.na(PeakNo_CC))
fungus.list <- separate(coculture.list, Sample_Ref, into = c("Ref_Culture", "Int_Culture"), sep = "v")
fungus <- unique(fungus.list$Ref_Culture)

for (i in length(fungus)) {
subset.list <- filter(full.list, Matched_con == paste0(fungus[i], "CON")) %>%
    select(Sample_Ref, PeakNo_con, PeakRatio, Metabolite_Effect) %>%
    mutate(logPeakRatio = log((PeakRatio/100) + 1))


heatmap_fvf <- ggplot(data = subset.list, aes(x=PeakNo_con, y=Sample_Ref)) + 
    geom_tile(aes(fill=logPeakRatio)) +
    scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                         midpoint = 0, limit = c(min(subset.list$logPeakRatio), 
                                                 max(subset.list$logPeakRatio)), 
                         space = "Lab", name="log(%PeakArea)") +
    labs(x="Peak Number", y="Coculture") +
    theme_grey(base_size=8)
heatmap_fvf
}

ggsave(heatmap_fvf,filename="heatmap_fvf.png",height=1.75,width=10.20,units="in",dpi=200)

#Testing

coculture.list <- filter(full.list, is.na(PeakNo_CC))
fungus.list <- separate(coculture.list, Sample_Ref, into = c("Ref_Culture", "Int_Culture"), sep = "v")
fungus <- unique(fungus.list$Ref_Culture)
fungus.of.interest <- paste0(substr(fungus[1],1,nchar(fungus[1])-3), "v")
    subset.list <- filter(full.list, grepl(fungus.of.interest, Sample_Ref)) %>%
        select(Sample_Ref, PeakNo_con, PeakRatio, Metabolite_Effect) %>%
        mutate(logPeakRatio = log((PeakRatio/100) + 1))
    
    
    heatmap_fvf <- ggplot(data = subset.list, aes(x=PeakNo_con, y=Sample_Ref)) + 
        geom_tile(aes(fill=logPeakRatio)) +
        scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                             midpoint = 0, limit = c(min(subset.list$logPeakRatio), 
                                                     max(subset.list$logPeakRatio)), 
                             space = "Lab", name="log(%PeakArea)") +
        labs(x="Peak Number", y="Coculture") +
        theme_grey(base_size=8)
    heatmap_fvf
    
    
    ## testing
    coculture.list <- filter(full.list, !is.na(PeakNo_CC)) %>%
        separate(Sample_Ref, into = c("Ref_Culture", "Int_Culture"), sep = "v")
    coculture.list$Ref_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", coculture.list$Ref_Culture, perl = TRUE)
    coculture.list$Int_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", coculture.list$Int_Culture, perl = TRUE)
    coculture.list$Matched_con <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", coculture.list$Matched_con, perl = TRUE)
     coculture.list <- arrange(coculture.list, Ref_Culture)
    fungus <- unique(coculture.list$Ref_Culture)

        subset.list <- filter(coculture.list, Ref_Culture == fungus[1], Matched_con == fungus[1]) %>%
            select(Int_Culture, PeakNo_con, Matched_con, PeakRatio) %>%
            mutate(logPeakRatio = log((PeakRatio/100) + 1.01))
        
        hm.title <- i
        heatmap_fvf <- ggplot(data = subset.list, aes(x=PeakNo_con, y=Int_Culture)) + 
            geom_tile(aes(fill=logPeakRatio)) +
            scale_fill_gradient2(low = "dark blue", high = "dark red", mid = "white", 
                                 midpoint = 0, limit = c(min(subset.list$logPeakRatio), 
                                                         max(subset.list$logPeakRatio)), 
                                 space = "Lab", name="log(%PeakArea)") +
            labs(x="Peak Number", y="Coculture") +
            theme_grey(base_size=8) +
            ggtitle(label = paste0("Interactions involving ", i))
        print(heatmap_fvf)
    
        
##Barchart testing
        
        p<- ggplot(Effect_Inductions, aes(x=Range, y=Proportion, fill=Category)) + 
            geom_bar(stat="identity", color="black", 
                     position=position_dodge())
        barchart <- p+labs(x="Elution Range (% MeCN)", y = "Proportion (%)")+
            theme_classic() +
            theme(axis.text.x=element_text(angle=45, hjust = 1))
        barchart
        ggsave(barchart,filename="barchart_fvf_tal.png",height=3.5,width=10.20,units="in",dpi=200)
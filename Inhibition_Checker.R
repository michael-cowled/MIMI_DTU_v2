#############################################
#14.Inhibition Checker: Verifies if one culture is inhibited
#############################################

##Note: removed due to bug, and not necessary for PHD

Inhibition_Checker <- function(df_Name, Inhibition_df, Coculture_Name) {
    
    #From an opened output file (might be possible just downstream of missing control peaks)
    match <- as.data.frame(table(df_Name$Matched_CON))
    names(match)[2] <- 'match'
    unmatch <- as.data.frame(table(df_Name$Sample_Ref))
    names(unmatch)[2] <- 'unmatch'
    combined <- merge(match, unmatch) %>%
        mutate(ratio = unmatch / match)
    
    #Then verify whether inhibition of control indicated
    
    if (nrow(combined) == 0) {
        combined <- match
        combined$unmatch[[1]] <- 0
        combined$ratio[[1]] <- 0
        combined$Inhibition[[1]] <- TRUE
    } else if (nrow(combined) == 1) {
        combined$Inhibition[[1]] <- TRUE
    }   else if (combined$match[1] <= 1 && combined$ratio[1] >= 3) {
        combined$Inhibition[[1]] <- TRUE
        combined$Inhibition[[2]] <- FALSE
    }   else if (combined$match[2] <= 1 && combined$ratio[2] >= 3) {
        combined$Inhibition[[2]] <- TRUE
        combined$Inhibition[[1]] <- FALSE
    }   else {
        combined$Inhibition[[2]] <- FALSE
        combined$Inhibition[[1]] <- FALSE
    }
    
    combined <- filter(combined, Inhibition == TRUE)
    if (nrow(combined) > 0) {
        combined$Coculture_Name <- Coculture_Name
        combined$Var1 <- as.character(combined$Var1)
        temp <- as.data.frame(c(combined[1,6], combined[1,5], combined[1,1]))
        rownames(temp) <- c("Coculture_Name", "Inhibition", "Inhibited_Culture")
        print(temp)
        Inhibition_df <- rbind(Inhibition_df, temp)
        print(Inhibition_df)
    }
    return(Inhibition_df)
    
}
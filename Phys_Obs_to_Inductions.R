induction.df$Ref_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", induction.df$Ref_Culture, perl = TRUE)
induction.df$Int_Culture <- gsub("(?<![0-9])([0-9])(?![0-9])", "0\\1", induction.df$Int_Culture, perl = TRUE)

induction.list <- merge(induction.df, Physical_Obs, by=c("Ref_Culture", "Int_Culture"), all.x=TRUE)

x <- filter(induction.list, Physical_Effect == 1, !is.na(induction.list$Num_Inductions))
sum(x$Num_Inductions)
mean(x$Num_Inductions)
qnorm(.95)*(sd(x$Num_Inductions)/sqrt(length(x)))

y <- filter(induction.list, Physical_Effect == 2, !is.na(induction.list$Num_Inductions))
sum(y$Num_Inductions)
mean(y$Num_Inductions)
qnorm(.95)*(sd(y$Num_Inductions)/sqrt(length(y)))

C <- filter(induction.list, Physical_Effect == "C", !is.na(induction.list$Num_Inductions))
sum(C$Num_Inductions)
mean(C$Num_Inductions)
qnorm(.95)*(sd(C$Num_Inductions)/sqrt(length(C)))

P <- filter(induction.list, Physical_Effect == "P", !is.na(induction.list$Num_Inductions))
sum(P$Num_Inductions)
mean(P$Num_Inductions)
qnorm(.95)*(sd(P$Num_Inductions)/sqrt(length(P)))

N <- filter(induction.list, Physical_Effect == "N", !is.na(induction.list$Num_Inductions))
sum(N$Num_Inductions)
mean(N$Num_Inductions)
qnorm(.95)*(sd(N$Num_Inductions)/sqrt(length(N)))

All <- filter(induction.list, !is.na(induction.list$Num_Inductions))

res <- t.test(P$Num_Inductions, All$Num_Inductions, alternative = "greater", var.equal = FALSE)
print(res$p.value)
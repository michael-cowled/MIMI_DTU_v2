# F22287:
F5vF15 <- read.csv("Testing Broad-Scale Interactions/OutputFiles/NT_FvF/F5vF15.CSV")

enh <- filter(F5vF15, Metabolite_Effect >3 & Metabolite_Effect <6, Sample_Ref != "F15CON", Matched_con != "F15CON")
supp <- filter(F5vF15, Metabolite_Effect <3, Matched_con != "F15CON")
mean(enh$RetTime_con)
mean(supp$RetTime_con)

res <- t.test(enh$RetTime_con, supp$RetTime_con, alternative = "less", var.equal = FALSE)
print(res$p.value)

# F22286:
F4vF15 <- read.csv("Testing Broad-Scale Interactions/OutputFiles/NT_FvF/F4vF15.CSV")

enh <- filter(F4vF15, Metabolite_Effect >3 & Metabolite_Effect <6, Sample_Ref != "F15CON", Matched_con != "F15CON")
supp <- filter(F4vF15, Metabolite_Effect <3, Matched_con != "F15CON")
mean(enh$RetTime_con)
mean(supp$RetTime_con)

res <- t.test(enh$RetTime_con, supp$RetTime_con, alternative = "two.sided", var.equal = FALSE)
print(res$p.value)

#induction between fvf and fva
ind_fvf <- filter(fvf_full.list, Metabolite_Effect == 6)
ind_fva <- filter(fva_full.list, Metabolite_Effect == 6)
res <- t.test(ind_fvf$RetTime_CC, ind_fva$RetTime_CC, alternative = "two.sided", var.equal = FALSE)
print(res$p.value)

#suppression between fvf and fva
supp_fvf <- filter(fvf_full.list, Metabolite_Effect == 2)
supp_fva <- filter(fva_full.list, Metabolite_Effect == 2)
res <- t.test(supp_fvf$RetTime_CC, supp_fva$RetTime_CC, alternative = "two.sided", var.equal = FALSE)
print(res$p.value)

#induction to supp. in fva
supp_fva <- filter(fva_full.list, Metabolite_Effect == 2)
ind_fva <- filter(fva_full.list, Metabolite_Effect == 6)
res <- t.test(supp_fva$RetTime_con, ind_fva$RetTime_CC, alternative = "two.sided", var.equal = FALSE)
print(res$p.value)

#induction to supp. in fvf
supp_fvf <- filter(fvf_full.list, Metabolite_Effect == 2)
ind_fvf <- filter(fvf_full.list, Metabolite_Effect == 6)
res <- t.test(ind_fvf$RetTime_CC, supp_fvf$RetTime_con, alternative = "two.sided", var.equal = FALSE)
print(res$p.value)

#enhancement to supp in fva
supp_fva <- filter(fva_full.list, Metabolite_Effect == 2)
enh_fva <- filter(fva_full.list, Metabolite_Effect == 4)
res <- t.test(enh_fva$RetTime_CC, supp_fva$RetTime_con, alternative = "two.sided", var.equal = FALSE)
print(res$p.value)

#enhancement to supp in fvf
supp_fvf <- filter(fvf_full.list, Metabolite_Effect == 2)
enh_fvf <- filter(fvf_full.list, Metabolite_Effect == 4)
res <- t.test(enh_fvf$RetTime_CC, supp_fvf$RetTime_con, alternative = "two.sided", var.equal = FALSE)
print(res$p.value)
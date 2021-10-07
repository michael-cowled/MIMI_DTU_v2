full_ret27 <- filter(full.list, RetTime_CC >=2 & RetTime_CC <=7 | RetTime_con >= 2 & RetTime_con <= 7)
full_ret24 <- filter(full.list, RetTime_CC >=2 & RetTime_CC <=4 | RetTime_con >= 2 & RetTime_con <= 4)
full_ret47 <- filter(full.list, RetTime_CC >=4 & RetTime_CC <=7 | RetTime_con >= 4 & RetTime_con <= 7)
enh <- filter(full.list, as.numeric(Metabolite_Effect) == 5 | as.numeric(Metabolite_Effect) == 4)
enh_27 <- filter(full_ret27, as.numeric(Metabolite_Effect) == 5 | as.numeric(Metabolite_Effect) == 4)
enh_24 <- filter(full_ret24, as.numeric(Metabolite_Effect) == 5 | as.numeric(Metabolite_Effect) == 4)
enh_47 <- filter(full_ret47, as.numeric(Metabolite_Effect) == 5 | as.numeric(Metabolite_Effect) == 4)
ind <- filter(full.list, as.numeric(Metabolite_Effect) > 5)
ind_27 <- filter(full_ret27, as.numeric(Metabolite_Effect) > 5)
ind_24 <- filter(full_ret24, as.numeric(Metabolite_Effect) > 5)
ind_47 <- filter(full_ret47, as.numeric(Metabolite_Effect) > 5)


length(full_ret27$Sample_Ref)/length(full.list$Sample_Ref)*100
length(full_ret24$Sample_Ref)/length(full.list$Sample_Ref)*100
length(full_ret47$Sample_Ref)/length(full.list$Sample_Ref)*100
length(enh_27$Sample_Ref)/length(enh$Sample_Ref)*100
length(enh_24$Sample_Ref)/length(enh$Sample_Ref)*100
length(enh_47$Sample_Ref)/length(enh$Sample_Ref)*100
length(ind_27$Sample_Ref)/length(ind$Sample_Ref)*100
length(ind_24$Sample_Ref)/length(ind$Sample_Ref)*100
length(ind_47$Sample_Ref)/length(ind$Sample_Ref)*100
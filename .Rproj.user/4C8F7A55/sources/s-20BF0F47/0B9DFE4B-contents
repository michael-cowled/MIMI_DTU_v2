X27v04$`UV Peaks` <- gsub("\\s*\\([^\\)]+\\)","",X27v04$`UV Peaks`)
X27v04 <- separate(X27v04, 'UV Peaks', c("UV1", "UV2", "UV3", "UV4", "UV5"), sep = "\r\n") %>%
  mutate(RT_h = RetTime + 0.5) %>%
  mutate(RT_l = RetTime - 0.5) %>%
  select(1:4, 11:15, 24:25)

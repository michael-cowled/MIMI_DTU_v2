#Before running functions, create a dataframe with the desried columns

df <- data.frame(PeakNo_CC = NA, RetTime_CC = NA, PeakArea_CC = NA,
                 PeakNo_CON = NA, RetTime_CON = NA, PeakArea_CON = NA,
                 UV_Count = NA, PeakRatio = NA)
print (df)

#Now to add the first lot of observations to row #1

df <- rbind(df, c(PeakNo_CC = CO27v04$Peak[i], RetTime_CC = CO27v04$RetTime[i], 
            PeakArea_CC = CO27v04$Area[i], PeakNo_CON = CON04$Peak[z], 
            RetTime_CON = CON04$RetTime[z], PeakArea_CON = CON04$Area[z], 
            UV_Count = uvcount, PeakRatio = ratio))

#Checking whether defining uvcount earlier works
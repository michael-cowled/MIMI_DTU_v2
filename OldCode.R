##Converting to RetTime columns to vectors for simplicity.

x <- as.vector(CO$RetTime)
y <- as.vector(df_04con$RetTime)

# x[n] from CO (the number to be tested)
# y from 04CON is the sequence to be tested

n <- length(x) #perhaps convert to nrow if performing this eventually on df.
i = 1
while (i < n+1) {
    z <- which(abs(y-x[i])==min(abs(y-x[i])))       #Finds the closest value to x[i] in y.
    if (x[i] < y[z] + 0.2 && x[i] > y[z] -0.2) {    #verifies that the closest value to x[i] in y is within +/- 0.2 min.
        cat("Peak#", i, x[i], "matches closest to Peak#", z, y[z], "in the control \n") 
        i <- i+1
    }   else  { 
        cat("Peak#", i, x[i], "is NOT present in the control \n") 
        i <- i+1
    }   
}
print("loop finished")


##Next to try without converting to vectors, i.e. as specific cells in a df.

# x = CO27v04
# y = CON04

n <- nrow(CO27v04) #converted to nrow from length
i = 1
while (i < n+1) {
    z <- which(abs(CON04$RetTime-CO27v04$RetTime[i])==min(abs(CON04$RetTime-CO27v04$RetTime[i])))
    if (CO27v04$RetTime[i] < CON04$RetTime[z] + 0.2 && CO27v04$RetTime[i] > CON04$RetTime[z] -0.2) {  
        cat("Peak#", i, "@", round(CO27v04$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON04$RetTime[z], digits =2), "min in the control \n") 
        i <- i+1
    }   else  { 
        cat("Peak#", i, "@", round(CO27v04$RetTime[i], digits =2), "min is NOT present in the control \n") 
        i <- i+1
    }   
}
print("loop finished")


##UV verification of 2 vectors: x/y (UV1:UV5 for 1 metabolite) [this works for missing values in both x and y]

uvcount = 0
j = 1
k = 1
while (k < 6) {
    j=1
    while (j < 6 && k < 6) {
        if (j==5)  {
            j <- j+1
            k <- k+1
            cat("UV count is", uvcount, "\n")
        }   else if (is.na(x[k])) {
            k <- k+1
        }   else if (is.na(y[j])) {
            j <- j+1
        }   else if (x[k] <= y[j] + 2 && x[k] >= y[j] - 2) {
            uvcount <- uvcount + 1
            j <- j+1
        }   else {
            j <- j+1
        }
    }
}
cat("Loop finished with uv count =", uvcount)

if (uvcount > 1) {
    print("metabolite confirmed")
}   else    {
    print ("not a match")
}

#Convert UV vector to Cells in a table.
#UVs are in columns 5-9, code rewritten to account for this
#Need to incorporate i and z from previous code to identify row numbers.
#Column numbers identified by j and k.
# x = CO27v04
# y = CON04

uvcount = 0
j = 5
k = 5 
while (k < 10) {
    j=5
    while (j < 10 && k < 10) {
        if (j==9)  {
            j <- j+1
            k <- k+1
            cat("UV count is", uvcount, "\n")
        }   else if (is.na(CO27v04[i,k])) {
            k <- k+1
        }   else if (is.na(CON04[z,j])) {
            j <- j+1
        }   else if (CO27v04[i,k] <= CON04[z,j] + 2 && CO27v04[i,k] >= CON04[z,j] - 2) {
            uvcount <- uvcount + 1
            j <- j+1
        }   else {
            j <- j+1
        }
    }
}
cat("Loop finished with uv count =", uvcount)

if (uvcount >= 1) {
    print("metabolite confirmed")
}   else    {
    print ("not a match")
}

#Now to combine the UV code into the Retention time matcher. Also added in the ratio of peak areas.

#But first, try encoding UV statements as a function:

UVcheck <- function() {
    uvcount <- 0
    j <- 5
    k <- 5 
    while (k < 10) {
        j<-5
        while (j < 10 && k < 10) {
            if (j==9)  {
                j <- j+1
                k <- k+1
            }   else if (is.na(CO27v04[i,k])) {
                k <- k+1
            }   else if (is.na(CON04[z,j])) {
                j <- j+1
            }   else if (CO27v04[i,k] <= CON04[z,j] + 2 && CO27v04[i,k] >= CON04[z,j] - 2) {
                uvcount <- uvcount + 1
                j <- j+1
            }   else {
                j <- j+1
            }
        }
    }
    if (uvcount >= 1) {
        cat("metabolite confirmed with uv count =", uvcount, "with a ratio of", round(ratio, digits =0), "% \n")
        return(uvcount)
    }   else    {
        cat("but UV does not match \n")
        return(uvcount)
    }   
}

#Now that the function works we will incorporate it into the RT checker

Coculture_df <- data.frame(PeakNo_CC = NA, RetTime_CC = NA, PeakArea_CC = NA,
                           PeakNo_CON = NA, RetTime_CON = NA, PeakArea_CON = NA,
                           UV_Count = NA, PeakRatio = NA)
n <- nrow(CO27v04)
i = 1
while (i < n+1) {
    z <- which(abs(CON04$RetTime-CO27v04$RetTime[i])==min(abs(CON04$RetTime-CO27v04$RetTime[i])))
    ratio = ((CO27v04[i,3]/CON04[z,3])*100) #Computes the ratio of peak areas as a %
    FinalCount <- UVcheck()
    if (CO27v04$RetTime[i] < CON04$RetTime[z] + 0.2 && CO27v04$RetTime[i] > CON04$RetTime[z] -0.2 && FinalCount > 0) {  
        cat("Peak#", i, "@", round(CO27v04$RetTime[i], digits =2), "min matches closest to Peak#", z, "@", round(CON04$RetTime[z], digits =2), "min in the control ")
        Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = CO27v04$Peak[i], RetTime_CC = CO27v04$RetTime[i], 
                                              PeakArea_CC = CO27v04$Area[i], PeakNo_CON = CON04$Peak[z], 
                                              RetTime_CON = CON04$RetTime[z], PeakArea_CON = CON04$Area[z], 
                                              UV_Count = FinalCount, PeakRatio = ratio))
        i <- i+1
    }   else  { 
        cat("Peak#", i, "@", round(CO27v04$RetTime[i], digits =2), "min is NOT present in the control \n")
        Coculture_df <- rbind(Coculture_df, c(PeakNo_CC = CO27v04$Peak[i], RetTime_CC = CO27v04$RetTime[i], 
                                              PeakArea_CC = CO27v04$Area[i], PeakNo_CON = NA, 
                                              RetTime_CON = NA, PeakArea_CON = NA, 
                                              UV_Count = NA, PeakRatio = NA))
        i <- i+1
    }   
}
print("loop finished")
write.csv(Coculture_df, paste0("C:/Users/Radicinol/Desktop/R/Testing Broad-Scale Interactions/OutputFiles/", COCulturename, ".CSV"), row.names = FALSE)

#Now to try the function on FP1927 and change to a function within a function and input samples to test

CON1 <- Read_Excel(CON1_Name)
CON2 <- Read_Excel(CON2_Name)
Coculture <- Read_Excel(Coculture_Name)


#Test UVCheck Function for MIMI

UVcheck <- function(CON1, Coculture, i, z) {
    uvcount <- 0
    j <- 5
    k <- 5 
    while (k < 10) {
        j<-5
        while (j < 10 && k < 10) {
            if (j==9)  {
                j <- j+1
                k <- k+1
            }   else if (is.na(Coculture[i,k])) {
                k <- k+1
            }   else if (is.na(CON1[z,j])) {
                j <- j+1
            }   else if (Coculture[i,k] <= CON1[z,j] + 2 && Coculture[i,k] >= CON1[z,j] - 2) {
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

#For CON2

UVCheck2 <- function(CON2, Coculture, i, z) {
    uvcount <- 0
    j <- 5
    k <- 5 
    while (k < 10) {
        j<-5
        while (j < 10 && k < 10) {
            if (j==9)  {
                j <- j+1
                k <- k+1
            }   else if (is.na(Coculture[i,k])) {
                k <- k+1
            }   else if (is.na(CON2[z,j])) {
                j <- j+1
            }   else if (Coculture[i,k] <= CON2[z,j] + 2 && Coculture[i,k] >= CON2[z,j] - 2) {
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

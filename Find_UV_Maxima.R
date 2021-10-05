raw <- read.csv("Simplified/Testing Broad-Scale Interactions/RawData_ava/8651CON_raw.CSV")
raw <- select(raw, -UV.Peaks)
uv <- read.csv("Simplified/Testing Broad-Scale Interactions/RawData_ava/8651CON_uv.CSV")


# Function for determining local maxima
find_peaks <- function (x, m = 20){ # 'm' is the stringency.
    shape <- diff(sign(diff(x, na.pad = FALSE)))
    pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
    })
    pks <- unlist(pks)
    pks
}

n <- ncol(uv)
for (i in 2:n) {
    print(i)
    a <- vector()
    b <- vector()
    Maxima <- find_peaks(uv[,i], m=20)
    print(Maxima)
    for (j in 1:length(Maxima)) {
        if (uv[Maxima[j], i] >= 5) {
            if (length(a) == 0) {
                a <- append(a, uv[Maxima[j], i])
                b <- paste0(b, uv[Maxima[j], 1], " (", round(uv[Maxima[j], i], 0), ")")
                print(b)
            }
            else if (uv[Maxima[j], i] >= 5 & abs(uv[Maxima[j], i] - as.numeric(a[length(a)])) > 1) {
                a <- append(a, uv[Maxima[j], i])
                b <- paste0(b, "\n", uv[Maxima[j], 1], " (", round(uv[Maxima[j], i], 0), ")")
                print(b) 
            }
        }
        
        if (length(b) >=1) {
            raw[i-1, 4] <- b
        }
        else {
            raw[i-1, 4] <- NA
        }
    }
}
names(raw)[4] <- "UV Peaks"


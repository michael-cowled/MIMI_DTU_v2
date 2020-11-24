UVSubtract1 <- function(CON1_UV, Coculture_UV, i, z) {
    k <- i + 1
    j <- z + 1
    SubtractedUV <- Coculture_UV[,k]-CON1_UV[,j]
    UV_Mean <- mean(SubtractedUV)
    print(UV_Mean)
    return(UV_Mean)
}

UVSubtract2 <- function(CON2_UV, Coculture_UV, i, z) {
    k <- i + 1
    j <- z + 1
    SubtractedUV <- Coculture_UV[,k]-CON2_UV[,j]
    UV_Mean <- mean(SubtractedUV)
    print(UV_Mean)
    return(UV_Mean)
}
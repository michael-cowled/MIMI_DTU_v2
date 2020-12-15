#This is a simple version of the effect categoriser with a single ratio input

Simple_Effect_Categoriser <- function(ratio) {
    if (ratio > -100 && ratio <= 20) {
        Effect <- 2
        return(Effect)
    }   else if (ratio > -20 && ratio < 20) {
        Effect <- 3
        return(Effect)
    }   else if (ratio >= 20 && ratio < 100) {
        Effect <- 4
        return(Effect)
    }   else if <- (ratio >= 100) {
        Effect <- 5
        return(Effect)
    }
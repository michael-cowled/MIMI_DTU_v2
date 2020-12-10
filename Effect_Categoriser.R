#The following piece of code adds a column that categorises the peak 
#areas into suppressions, and enhancements

Effect_Categoriser <- function(Refined_Coculture_df, Coculture_Name) {

    RowNo <- 1
    Rows <- nrow(Refined_Coculture_df)
    
    #A df is created to list the effects corresponding to matched peaks
    
    Metabolite_effect_df <- setNames(data.frame(matrix(ncol=1, nrow=Rows)), 
                                     c("Metabolite_Effect"))
    
    while (RowNo <= Rows) {
        if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
            is.na(Refined_Coculture_df[RowNo, 11])) {
            Metabolite_effect_df[RowNo, 1] <- 6 #Induction
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] > -100
                     && Refined_Coculture_df[RowNo, 11] <= -20) {
            Metabolite_effect_df[RowNo, 1] <- 2 #Suppression
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] > -20
                     && Refined_Coculture_df[RowNo, 11] < 20) {
            Metabolite_effect_df[RowNo, 1] <- 3 #Little to No Change
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] >= 20
                     && Refined_Coculture_df[RowNo, 11] < 100) {
            Metabolite_effect_df[RowNo, 1] <- 4 #Enhancement
            RowNo <- RowNo + 1
        }   else if (Refined_Coculture_df[RowNo, 1] == Coculture_Name && 
                     Refined_Coculture_df[RowNo, 11] >= 100) {
            Metabolite_effect_df[RowNo, 1] <- 5 #Complete Suppression
            RowNo <- RowNo + 1
        }
    }
    Refined_Coculture_df <- 
        cbind(Refined_Coculture_df, Metabolite_effect_df)
    return(Refined_Coculture_df)
}
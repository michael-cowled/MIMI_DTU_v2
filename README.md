# Microbial_Interactions    ## [1] "Date is: Tue 07 Sep 2021"

## R Packages required to be loaded in:

    library(dplyr)
    library(tidyr)
    library(readxl)

## NovaC files

NovaC is an unpublished data extraction and formatting tool for Agilent
Chemstation LCMS/HPLC files owned by Macquarie University and created by
A. Prof. Andrew Piggott (Molecular Sciences/Faculty of Science &
Engineering).

The raw data was extracted from Agilent Chemstation LCMS/HPLC using the
macro “NovaC.mac” (Macquarie Unviersity) which generates an Excel file
with the following name: "002-D1F-D7-\_F7vF2.D\_NovaC", as an example.

See example NovaC files in folder titled: “Example\_NovaC\_Data”.

File names were manually tidied up to remove unnecessary prefixes or
suffixes: “F7vF2”.

This name tidying and the creation of the Interaction\_Matrix below, are
the only manual steps performed.

## The Interaction Matrix

A table in excel format with the heading names as follows is required to
be loaded in by the user. The matrix should contain all interactions to
be investigated.

    Interaction_Matrix <- read_excel("Interaction_Matrix.xlsx")
    head(Interaction_Matrix)

    ## # A tibble: 6 x 3
    ##   CON1  CON2  Coculture
    ##   <chr> <chr> <chr>    
    ## 1 F1CON F2CON F1vF2    
    ## 2 F1CON F3CON F1vF3    
    ## 3 F1CON F4CON F1vF4    
    ## 4 F1CON F5CON F1vF5    
    ## 5 F1CON F6CON F1vF6    
    ## 6 F1CON F7CON F1vF7

## MIMI - Script \#1/2

The following functions are to be pre-loaded prior to use of the main
function, MIMI()

1.  **ReadExcel:** Reads in the NovaC files corresponding to a
    particular row number in the Interaction\_Matrix.
2.  **CheckUVCount:** Compares the top 5 UV maxima of a matched peak in
    the control of interest and the coculture.
3.  **ReadUV:** Reads in the UV spectral data (abs vs. wavelength) for
    the corresponding 3 samples being compared.
4.  **SubtractUV:** Subtracts the UV spectrum of the control of interest
    from the coculture, and takes the mean(Abs).
5.  **RowBinder:** Binds the matched peak found in PeakMatcher to a df
    named cc.df.
6.  **CalcRatio:** Calculates the %Enhancement/Suppresion compared to
    control levels.
7.  **PeakMatcher:** Matches and verifies peaks from the control of
    interest to peaks in the coculture, utilising a combination of
    retention time, number of matching UV maxima (UVcheck) and the means
    of the subtracted UV spectra (UVsubtract).
8.  **ConConsolidator:** Consolidates the outcome table to match peaks
    from coculture to a single peak in a control. Compares subtracted UV
    spectra to make decisions based on double matching.
9.  **EffectCategoriser:** Characterises the Peak Area ratio as an
    effect to the metabolite in the coculture (induction, suppression,
    etc.).

## MIMI - Script \#2/2

The following functions are to be pre-loaded prior to use of the main
function, MIMI2()

1.  **SimpleEffectCharacteriser:** Uses the principles of the
    EffectCategoriser function to categorise based on a single PeakRatio
    input.
2.  **DoublePeakChecker:** Finds and reports instances of double peak
    matching.
3.  **IdentifyBadPeak:** Decides on which of the doubley assigned peaks
    to remove.
4.  **RemoveDoubleyAssignedPeaks:** Checks for a peak in the control
    being assigned to more than one peak in the coculture. Uses the
    number of matching UV maxima and/or the subtracted UV spectra to
    make decisions as to which peak is a better match.
5.  **PeakAssigner:** Adds the assignment of the matched non UV peak
    from MatchNonUVs.
6.  **MatchNonUVs:** Tentatively assigns matched peaks as non UVs (or as
    distorted UVs) if matching the conditions.
7.  **RowBinder2:** A variant of RowBinder used for adding missing
    control peaks.
8.  **FindMissingcontrolPeaks:** Adds in unassigned peaks from the
    controls to provide a single, unified table.

Please see the full script to load in the whole suite of functions and
review annotations for extra details. The final main functions,
**MIMI()** and **MIMI2()**, combines the use of all 17 functions in
providing the final tidied output files.

## The output file

An example output file is as follows:

    ##   Sample_Ref PeakNo_CC RetTime_CC PeakArea_CC   PercArea Matched_con PeakNo_con RetTime_con PeakArea_con UV_Count
    ## 1      F1vF3         1   0.574947   532.11786 18.0381983       F1CON          1    0.584723    265.30689        1
    ## 2      F1vF3         2   0.823738   144.50023  4.8983956       F1CON          2    0.850788     78.56687        2
    ## 3      F1vF3         3   1.039050   943.08105 31.9693895       F1CON          3    1.060212    536.97082        4
    ## 4      F1vF3         4   1.943127   353.14883 11.9713492       F1CON          4    1.970453     33.35677        2
    ## 5      F1vF3         5   2.600197    94.46783  3.2023533       F1CON          5    2.609808     30.68041        2
    ## 6      F1vF3         6   2.714563    18.00540  0.6103629       F1CON          6    2.723022     18.68216        1
    ##   Subtracted_UV_Mean  PeakRatio Metabolite_Effect
    ## 1         0.13295245 100.566924                 5
    ## 2        -0.43382285  83.920048                 4
    ## 3        -0.03814868  75.629850                 4
    ## 4        -0.22120232 958.702255                 5
    ## 5         0.05347423 207.909232                 5
    ## 6        -0.36560138  -3.622473                 3

The resulting output file contains the following columns containing:

1.  **Sample\_Ref:** The reference sample - the default is the
    coculture, unless the peak from a control is not present.
2.  **PeakNo\_CC:** The peak number as present in the coculture sample.
3.  **RetTime\_CC:** The retention time (in min) for the peak in the
    coculture sample.
4.  **PeakArea\_CC:** The peak area for the peak in the coculture
    sample.
5.  **PercArea:** The percentage the peak area associated with this peak
    accounts for compared to the rest of the peaks in the sample.
6.  **Matched\_con:** The control sample with a matching peak to the
    sample reference.
7.  **PeakNo\_con:** The peak number as present in the matched control
    or reference control sample.
8.  **RetTime\_con:** The retention time (in min) for the peak in the
    control sample.
9.  **PeakArea\_con:** The peak area for the peak in the control sample.
10. **UV\_Count:** The number of UV maxima that match between the
    control and coculture peaks.
11. **Subtracted\_UV\_Mean:** The mean of the subtracted UV spectra (in
    Absorbance units) of the control from the coculture.
12. **PeakRatio:** The relative ratio (in %) of the reference sample
    compared to the control.
13. **Metabolite\_Effect:** A categorised effect on the metabolite in
    the reference sample as determined by the PeakRatio. These effects
    are: *Complete Suppression:* Peak Area = -100% | *Suppression:*
    -100% &lt; Peak Area &lt;= -20% | *Little to No Change:* -20% &lt;
    Peak Area &lt; 20% *Enhancement:* 20 &lt;= Peak Area &lt; 100% |
    *Major enhancement:* Peak area &gt;=100% | *Induction:* NA (No peak
    area to compare to in CON)

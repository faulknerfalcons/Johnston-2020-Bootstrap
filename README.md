# Johnston-2020-NewPhyt
## A bootstrap approach is a superior statistical method for the comparison of cell-to-cell movement data
Plasmodesmata are an increasing focus of plant research, and plant physiologists frequently aim to understand the dynamics of intercellular movement and plasmodesmal function. For this, experiments that measure the spread of GFP between cells are commonly performed to indicate whether plasmodesmata are more open or closed in different conditions or in different genotypes. 
We propose that cell-to-cell movement data sets are better analysed by a bootstrap method that tests the null hypothesis that means (or medians) are the same between two conditions, instead of a Mann-Whitney-Wilcoxon test.  We found that that with hypothetical distributions similar to cell-to-cell movement data, the Mann-Whitney-Wilcoxon produces a false positive rate of 17% while the bootstrap method maintains a false positive at the set rate of 5% under the same circumstances.  Here we present this finding, as well as our rationale, an explanation of the bootstrap method and an R script for easy use. We have further demonstrated its use on published datasets from independent laboratories.

### RawData.xlsx
The data in Fig. 1 comes from (b) Figure S2 of Cheval et al. (2020), (c) Figure 2d Diao et al. (2018), (d) Figure 2c Diao et al. (2018) under use of the CC BY 4.0 licence.

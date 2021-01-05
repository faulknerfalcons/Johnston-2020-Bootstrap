## A bootstrap approach is a superior statistical method for the comparison of non-normal data with differing variances
#### New Phytologist (2020) doi: [10.1111/nph.17159](https://doi/org/10.1111/nph.17159)

#### Quick Use
`source("https://raw.githubusercontent.com/faulknerfalcons/Johnston-2020-Bootstrap/master/medianBootstrapToolbox.R")`

`medianBootstrap(data1, data2)`

#### Abstract
Phytologists rely on experimentally perturbing plants and monitoring the responses. Frequentist statistics are used to ascertain the probability that an observed difference between conditions was due to chance (a p value). When data are not normal and have differing variances, we propose data sets are better analysed by a bootstrap method that tests the null hypothesis that means (or medians) are the same between two conditions, instead of the commonly used  Mann-Whitney-Wilcoxon test. 

We illustrate this with data from the cell-to-cell movement of GFP through plasmodesmata. We found that that with hypothetical distributions similar to cell-to-cell movement data, the Mann-Whitney-Wilcoxon produces a false positive rate of 17% while the bootstrap method maintains a false positive at the set rate of 5% under the same circumstances.  Here we present this finding, as well as our rationale, an explanation of the bootstrap method and an R script for easy use. We have further demonstrated its use on published datasets from independent laboratories.

This GitHub repository is archived at Zenodo

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4326968.svg)](https://doi.org/10.5281/zenodo.4326968)

### RawData.xlsx
The data in Fig. 1 comes from (b) Figure S2 of Cheval et al. (2020), (c) Figure 2d Diao et al. (2018), (d) Figure 2c Diao et al. (2018) under use of the CC BY 4.0 licence.

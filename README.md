Sufficient Dimension Reduction for gene expression data.

Sufficient dimension reduction is a supervised learnig method for reducing the number of random vairables into a set of principal variables along with the notion of sufficiency. 

Before applyig SDR, we first need to do a variable selection of genes using ANOVA testing. new_data_anova.rda is the data used in the program. The data is originally acquired from R data package 'rcellminerdata' from Bioconductor (https://www.bioconductor.org/packages/release/data/experiment/html/rcellminerData.html). I have used copy number, micro RNA, mRNA and proteomics. 

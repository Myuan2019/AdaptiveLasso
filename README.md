# AdaptiveLasso

There are two files provided for LASSO_FBAT analysis: 

The "LASSO_FBAT sample.Rdata" contains a set of simulated example of family-structured data. There are two variables contained in the Rdata file.

Variable "sample" is a three-dimention array with the first dimention to be the gene index, the second dimention to be the family index.
The third dimention of variable "sample" is  respectively: gene score of parent 1, gene score of parent 2, gene score of the offspring 
and expectation of gene score based on parental information (i.e. E(X|P)).

Variable "diserisk" is a two-dimention array contain the disease status of the parents.


The "code for LASSO_FBAT analysis.R" contains the R code for analysis.

Please load the sample Rdata file first and then do the analysis. 

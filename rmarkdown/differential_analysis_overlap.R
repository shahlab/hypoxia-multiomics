###
###In the following we assume we already have the results objects from the differential analyses of two experiments. 
###One object is res_exp1 and the other is res_exp2. We assume that we have well-calibrated p-value distributions, 
###and the results objects are in the format provided by deseq2
###

##First, to visually get a feel for the overlap between the differential hits, look at the conditional p-value distributions
hist(res_exp2$pvalue[which(rownames(res_exp2) %in% rownames(res_exp1)[which(res_exp1$padj < 0.1)])], 
     breaks = 40, freq = FALSE)
hist(res_exp2$pvalue[-which(rownames(res_exp2) %in% rownames(res_exp1)[which(res_exp1$padj < 0.1)])], 
     breaks = 40, freq = FALSE)


##now use the qvalue package to estimate pi_0 and label genes as shared differential at the desired fdr level
library(qvalue)

qobj <- qvalue(p = res_exp2$pvalue[which(rownames(res_exp2) %in% rownames(res_exp1)[which(res_exp1$padj < 0.1)])], 
                                        fdr.level = 0.1, pi0.method = "bootstrap")

#note: it is possible for the estimated pi_0 to be discordant (bigger or smaller) with what we intuitively expect 
#based on the p-val distributions. In this case, look at the pi_0 values produced for different values of the lambda parameter

shared_differential_genes <- res_exp2[which(rownames(res_exp2) %in% 
                                             rownames(res_exp1)[which(res_exp1$padj < 0.1)])[
                                                          which(qobj$significant == TRUE)], ]

##get a p-value for the overlap via permutations
null_distribution <- replicate(10000, {
  pi0_random <- pi0est(p = sample(res_exp2$pvalue, 
                                  length(which(rownames(res_exp2) %in% rownames(res_exp1)[which(res_exp1$padj < 0.1)]))), 
                       pi0.method = "bootstrap")$pi0
  1- pi0_random
})






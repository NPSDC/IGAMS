# IGAMS - Integrated Models for Gene Expression and DNA Methylation for Cancer Stage Prediction

This repository contains the implementation of the models demonstrated in the paper
"**Integrative analysis of DNA methylation and gene expression in Papillary Renal Cell Carcinoma**".

We have tried two different MKL (Multiple Kernel Learning) approaches - Group Lasso and BEMKL for integrating data across different platforms, the code for which has been taken and adapted from the following two repositories:

Bayesian Efficient Multiple Kernel Learning - <https://github.com/mehmetgonen/bemkl>  
Group Lasso - <https://github.com/mehmetgonen/gsbc>


Always stay in the home directory in order to run any script i.e R environmnent should be set in home directory.
meth_train.R
rna_meth_train.R
rna_train.R   
These files contain the examples of how to make use of the pipeline and the steps that we followed as a part of our analysis.
Note we expect the data from the different platforms as a data.frame as input, everything else would be taken care by the pipeline.

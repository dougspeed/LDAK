LDAK

July 2023. LDAK now includes TetraHer and QuantHer, methods for estimating family heritability for binary and quantitative traits. We will shortly upload a preprint providing methodological details.

June 2022. We have added LDAK-GBAT, our new tool for gene-based association testing using GWAS summary statistics. Details of the tool can be found in our our paper LDAK-GBAT: fast and powerful gene-based association testing using summary statistics, published in the American Journal of Human Genetics, with scripts from the paper provided in Takiy's GitHub.

January 2022. We have just released version 5.2 of LDAK. This version includes Quick PRS, a super-fast way to construct state-of-the-art SNP-based prediction models for complex human traits that requires only summary statistics. You can obtain LDAK v5.2 from Downloads.

July 2021. Our paper, Improved genetic prediction of complex traits from individual-level data or summary statistics, has now been published in Nature Communications. The LDAK software has been updated to include the seven new Prediction tools described in the paper.

March 2020. Our paper Evaluating and improving heritability models using summary statistics has been published in Nature Genetics (for a free version click here). This work represents a major update to SumHer, our tool for performing heritability analyses using summary statistics. In particular, it provides methods for comparing heritability models and for estimating the selection-related parameter alpha.
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

LDAK is a software package for analysing association study data. In total, it has over 30 different functions (see Main Arguments for a full list). However, most likely, you will be interested in one of the following:

Testing predictors for association with a phenotype (either individually or jointly).

Estimating SNP heritability by analysing individual-level data

Estimating SNP heritability, heritability enrichments, genetic correlation and the selection-related parameter alpha by analysing summary statistics

Constructing SNP-based prediction models (polygenic risk scores)
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _

Central to LDAK is the Heritability Model, which describes the assumed distribution of heritability across the genome. In human genetics, it is common to assume that each SNP contributes equally. However, we have shown that it is better to use models where heritability depends on factors such as allele frequency and linkage disequilibrium (see Publications for more details).

If you are new to LDAK, I suggest you first obtain the latest version from Downloads, then read the Advice. If you have any problems, or spot any errors on this website, please contact me (see Help for details).

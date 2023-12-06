# Joint-Lassosum 

Welcome! Here we have source code for PolyGenic risk Score method Joint-Lassosum.

## Overview

Joint-Lassosum is a statistical method designed to enhance the portability of polygenic scores (PGS) across different ancestral groups. Polygenic scores are numerical indicators used to predict an individual's susceptibility to a particular disease. However, PGS developed in one population may not be as accurate when applied to individuals from diverse ancestral backgrounds. Joint-Lassosum addresses this issue by incorporating data from two distinct populations, resulting in a more universally applicable PGS, which can especially benefit the under-represented population.

For more information on the method and its implementation, please refer to the associated paper.

## Citation

If you use the Joint-Lassosum method or code for your research, please cite our paper "Evaluating and Improving Health Equity and Fairness of Polygenic Scores"

## Basic Usage

New users may start with [`fit_JLS.R`](/JLS_basic/code/fit_JLS.R) - Given GWAS information and reference population genotype (for calculating gene-gene correlation), it fits Joint-Lassosum. Please make sure your data is in the same format as our [toy example](/JLS_basic/data/).

After fitting the JLS models, essentially calculating the high-dimensional Lasso-type regression coefficients, we can use [`evaluate_PGS.R`](/JLS_basic/code/evaluate_PGS.R) to predict the PGS of individuals given their genotype. We present the data management steps for the ease of use.

## About Reference Panel/Population

(Joint-)Lassosum-type algorithms leverage the SNP correlation information of different populations. In addition to the GWAS summary statistics, JLS also requires the users to supply a set of plink files (containing individual genotype information) to calculate the correlation between SNPs within each LD block. In our toy example, they are the files starting with "CEU-" and "YRI-" in the [data folder](/JLS_basic/data/). Note we split the whole genome information into chromosome chunks to facilitate parallelization. 

The individuals do not need to be those involved in the GWAS study---we even allow them to be simulated synthetic "individuals" so long as the SNP information is preserved. We suggest the users download the original data from 1000 Genomes Project in their cluster computer and sequentially run our [simualtion pipeline](/Generate_reference_populaton/) to generate a "synthetic" individual in the required format. 








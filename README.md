# Joint-Lassosum 

Source code for Polygenic Risk Score methods Joint-Lassosum.

## Overview

Joint-Lassosum is a statistical method designed to enhance the portability of polygenic scores (PGS) across different ancestral groups. Polygenic scores are numerical indicators used to predict an individual's susceptibility to a particular disease. However, PGS developed in one population may not be as accurate when applied to individuals from diverse ancestral backgrounds. Joint-Lassosum addresses this issue by incorporating data from two distinct populations, resulting in a more universally applicable PGS, which can especially benefit the under-represented population.

For more information on the method and its implementation, please refer to the associated paper.

## Citation

If you use the Joint-Lassosum method or code for your research, please cite our paper "Evaluating and Improving Health Equity and Fairness of Polygenic Scores"

## Usage

New users may start with [`fit_JLS.R`]/JLS_basic/code/fit_JLS.R) - Given GWAS information and reference population genotype (for calculating gene-gene correlation), it fits Joint-Lassosum. Please make sure your data is in the same format as those in (https://github.com/terrytianyuzhang/JointLassosum/tree/main/JLS_basic/data)










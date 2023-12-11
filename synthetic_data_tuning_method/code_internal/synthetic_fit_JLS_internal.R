gc()

library(data.table)
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(R.utils)
# library(pROC) # for AUC 
# library(pryr) # check memory useage

source("/raid6/Tianyu/PRS/sharable_synthetic_tuning/code/generate_synthetic_data_function.R")
source('/raid6/Tianyu/PRS/sharable/code/JLS_function.R')
source('/raid6/Tianyu/PRS/sharable/code/JLS_function_supp.R')
Rcpp::sourceCpp("/raid6/Tianyu/PRS/sharable/code/JLS_function.cpp")

### NEED TO TELL ME THE SYNTHETIC HYPERPARAMETER SO THAT I CAN FIND THE FILES
synthetic_population_weight_one <- 0.5 ###THIS NEED TO BE ONE OF THE CANDIDATE "GAMMA"
synthetic_l1_penalty_one <- 0.02236068###THIS NEED TO BE ONE OF THE CANDIDATE "LAMBDA"

para_tuning_result_folder <- '/raid6/Tianyu/PRS/sharable_synthetic_tuning/result_internal/' ###PLACE TO STORE PARAMETER TUNING RESULTS
large_population_type <- 'CEU' #or EUR, for determining LD block boundaries
small_population_type <- 'YRI' #or AFR

###HYPERPARAMETER CANDIDATES
JLS_population_weight <- c(0.2, 0.5, 0.8) #gamma parameter in the paper
# JLS_l1_penalty <- exp(seq(log(0.007), log(0.05), length.out=5))
JLS_l1_penalty <- exp(seq(log(0.01), log(0.05), length.out = 5)) ##lambda
JLS_shrinkage <- c(0.9) #s
chromosome <- 1:22###usually you need to change this to 1:22 

###LOCATION OF THE SYNTHETIC FIT RESULTS
synthetic_JLS_result_folder <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                         'synthetic_refit_',
                                                                         synthetic_population_weight_one,
                                                                         synthetic_l1_penalty_one,
                                                                         '/')
print(synthetic_JLS_result_folder)


### FILE LOCATION, THIS MATCHES HOW THE SYNTHETIC DATA SETS ARE GENERATED AND DO NOT NEED TO MODIFIED
JLS_result_prefix <- paste0(synthetic_JLS_result_folder, 'synthetic_JLS_result_weight_is')
large_population_GWAS_folder <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                          paste0(large_population_type, '_synthetic_GWAS_'),
                                                                          synthetic_population_weight_one,
                                                                          synthetic_l1_penalty_one,
                                                                          '/')
large_population_GWAS_file <- paste0(large_population_GWAS_folder, 'gwas_all.PHENO1.glm.logistic.hybrid')
small_population_GWAS_folder <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                          paste0(small_population_type, '_synthetic_GWAS_'),
                                                                          synthetic_population_weight_one,
                                                                          synthetic_l1_penalty_one,
                                                                          '/')
small_population_GWAS_file <- paste0(small_population_GWAS_folder, 'gwas_all.PHENO1.glm.logistic.hybrid')
large_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/CEU/CHR/CEU-chr'
small_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/YRI/CHR/YRI-chr'


###OTHER METADATA

# small_population_type <- 'ASN'


###GIVEN THE ABOVE INFORMATION, FIT THE MODEL FOR ONE TIME
# JLS_population_weight_one <- JLS_population_weight[1]

#####BASIC USAGE FOR ONE GIVEN POPULATION WEIGHT PARAMETER, THE OTHER TWO PARAMETERS CAN BE VECTORS
# JLS_train(JLS_population_weight_one = JLS_population_weight_one,
#           JLS_l1_penalty = JLS_l1_penalty,
#           JLS_shrinkage = JLS_shrinkage,
#           large_population_GWAS_file = large_population_GWAS_file,
#           small_population_GWAS_file = small_population_GWAS_file,
#           large_population_reference_prefix = large_population_reference_prefix,
#           small_population_reference_prefix = small_population_reference_prefix,
#           JLS_result_prefix = JLS_result_prefix,
#           large_population_type = large_population_type,
#           small_population_type = small_population_type,
#           chromosome = chromosome)

##PARALLEL OVER DIFFERENT WEIGHT PARAMETER GAMMA
mclapply(JLS_population_weight,
         FUN = JLS_train,
         JLS_l1_penalty = JLS_l1_penalty,
         JLS_shrinkage = JLS_shrinkage,
         large_population_GWAS_file = large_population_GWAS_file,
         small_population_GWAS_file = small_population_GWAS_file,
         large_population_reference_prefix = large_population_reference_prefix,
         small_population_reference_prefix = small_population_reference_prefix,
         JLS_result_prefix = JLS_result_prefix,
         large_population_type = large_population_type,
         small_population_type = small_population_type,
         chromosome = chromosome,
         num_parallel_chr = 8, ###change this to a larger number for palatalization over chromosomes.
         mc.cores = 2, ###change this to a larger number for palatalization over weight parameters.
         mc.preschedule = F, 
         mc.silent = F)

####NOW THE FITTING RESULTS ARE IN THE SPECIFIED LOCATION

# large_population_GWAS_two_chr <- large_population_GWAS[large_population_GWAS$`#CHROM` %in% c(21,22), ]
# small_population_GWAS_two_chr <- small_population_GWAS[small_population_GWAS$`#CHROM` %in% c(21,22), ]
# large_population_GWAS_two_chr_file <- '/raid6/Tianyu/PRS/sharable/large_population_GWAS_two_chr'
# small_population_GWAS_two_chr_file <- '/raid6/Tianyu/PRS/sharable/small_population_GWAS_two_chr'
# fwrite(large_population_GWAS_two_chr, file = large_population_GWAS_two_chr_file)
# fwrite(small_population_GWAS_two_chr, file = small_population_GWAS_two_chr_file)


#

# source("Lassosum_function/parseselect.R")
# source("Lassosum_function/parseblocks.R")
# source("Lassosum_function/ncol.bfile.R")
# source("Lassosum_function/nrow.bfile.R")
# source("Lassosum_function/read.table2.R")
# source("Lassosum_function/selectregion.R")
# source("Lassosum_function/parse.pheno.covar.R")
# source("Lassosum_function/splitgenome.R")
# source("Lassosum_function/validation.R")
# source("Lassosum_function/merge.mylassosum.R")

# large_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid'
# small_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/YRI.TRN.PHENO1.glm.logistic.hybrid'
# large_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/CEU/CHR/CEU-chr'
# small_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/YRI/CHR/YRI-chr'


# large_population_GWAS_file <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/large_population_GWAS_two_chr'
# small_population_GWAS_file <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/small_population_GWAS_two_chr'
# large_population_reference_prefix <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/CEU-chr'
# small_population_reference_prefix <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/YRI-chr'
# gene_ID_chromosome_file <- '/raid6/Ron/prs/data/bert_sample/CEU.TRN/CEU.TRN.pvar'


###########generate boostrap samples and split them into training and validation######
gc()

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(R.utils)
library(snpStats)
library(wordspace)
library(foreach)

#####
source("/raid6/Tianyu/PRS/sharable_synthetic_tuning/code/generate_synthetic_data_function.R")

JLS_population_weight_one <- 0.5 ###THIS NEED TO BE ONE OF THE CANDIDATE "GAMMA"
JLS_l1_penalty_one <- 0.02236068###THIS NEED TO BE ONE OF THE CANDIDATE "LAMBDA"

##the result of one JLS fit is the input of this pipeline
JLS_result_folder <- '/raid6/Tianyu/PRS/sharable/result_internal/' ###POINT THE CODE TO THE OUTPUT OF JLS FITTING RESULT. AFTER RUNNING fit_JLS.R
large_population_type <- 'CEU'
small_population_type <- 'YRI'

synthetic_large_population_prefix <- '/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/CEU-20K/CEU-20K'
synthetic_small_population_prefix <- '/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/YRI-4K/YRI-4K'
synthetic_large_population_prefix_by_chr <- '/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/CEU-20K/CHR/CEU-20K-chr'
synthetic_small_population_prefix_by_chr <- '/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/YRI-4K/CHR/YRI-4K-chr'

# large_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/CEU/CHR/CEU-chr'
# small_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/YRI/CHR/YRI-chr'

large_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid'
small_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/YRI.TRN.PHENO1.glm.logistic.hybrid'

large_population_GWAS_case_proportion <- 0.5
small_population_GWAS_case_proportion <- 0.5
chrs <- 1:22



###OUTPUT FOLDER
para_tuning_result_folder <- '/raid6/Tianyu/PRS/sharable_synthetic_tuning/result_internal/' ###PLACE TO STORE PARAMETER TUNING RESULTS

#####STEP 1: DETERMINE THE PGS FOR EACH INDIVIDUAL IN THE SYNTHETIC POPULATION
predict_PGS_given_coefficient(JLS_result_folder = JLS_result_folder,
                              para_tuning_result_folder = para_tuning_result_folder,
                              population_type = large_population_type,
                              synthetic_population_prefix = synthetic_large_population_prefix,
                              JLS_population_weight_one = JLS_population_weight_one,
                              JLS_l1_penalty_one = JLS_l1_penalty_one)

predict_PGS_given_coefficient(JLS_result_folder = JLS_result_folder,
                              para_tuning_result_folder = para_tuning_result_folder,
                              population_type = small_population_type,
                              synthetic_population_prefix = synthetic_small_population_prefix,
                              JLS_population_weight_one = JLS_population_weight_one,
                              JLS_l1_penalty_one = JLS_l1_penalty_one)

#####STEP 2: ASSIGN LABEL PROVIDED PGS
synthetic_label_given_PGS(JLS_result_folder = JLS_result_folder,
                          GWAS_file = small_population_GWAS_file,
                          population_type = small_population_type,
                          para_tuning_result_folder = para_tuning_result_folder,
                          case_proportion = small_population_GWAS_case_proportion,
                          JLS_population_weight_one = JLS_population_weight_one,
                          JLS_l1_penalty_one = JLS_l1_penalty_one)

synthetic_label_given_PGS(JLS_result_folder = JLS_result_folder,
                          GWAS_file = large_population_GWAS_file,
                          population_type = large_population_type,
                          para_tuning_result_folder = para_tuning_result_folder,
                          case_proportion = large_population_GWAS_case_proportion,
                          JLS_population_weight_one = JLS_population_weight_one,
                          JLS_l1_penalty_one = JLS_l1_penalty_one)
set.seed(2019)
# num_fold <- 5
# chrs <- 21:22
# population_type <- 'CEU'
# GWAS_file <- small_population_GWAS_file
# synthetic_population_prefix_by_chr<- synthetic_small_population_prefix_by_chr

mclapply(chrs, split_train_validation, 
         population_type = large_population_type,
         para_tuning_result_folder = para_tuning_result_folder,
         synthetic_population_prefix_by_chr = synthetic_large_population_prefix_by_chr,
         GWAS_file = large_population_GWAS_file, 
         JLS_population_weight_one = JLS_population_weight_one,
         JLS_l1_penalty_one = JLS_l1_penalty_one,
         mc.cores = 4, mc.preschedule = FALSE)
mclapply(chrs, split_train_validation, 
         population_type = small_population_type,
         para_tuning_result_folder = para_tuning_result_folder,
         synthetic_population_prefix_by_chr = synthetic_small_population_prefix_by_chr,
         GWAS_file = small_population_GWAS_file, 
         JLS_population_weight_one = JLS_population_weight_one,
         JLS_l1_penalty_one = JLS_l1_penalty_one,
         mc.cores = 4, mc.preschedule = FALSE)

calculate_synthetic_GWAS(population_type = large_population_type,
                         para_tuning_result_folder = para_tuning_result_folder,
                         synthetic_population_prefix_by_chr = synthetic_large_population_prefix_by_chr,
                         GWAS_file = large_population_GWAS_file,
                         JLS_population_weight_one = JLS_population_weight_one,
                         JLS_l1_penalty_one = JLS_l1_penalty_one,
                         chrs = chrs,
                         num_chr_parallel = 4)
calculate_synthetic_GWAS(population_type = small_population_type,
                         para_tuning_result_folder = para_tuning_result_folder,
                         synthetic_population_prefix_by_chr = synthetic_small_population_prefix_by_chr,
                         GWAS_file = small_population_GWAS_file,
                         JLS_population_weight_one = JLS_population_weight_one,
                         JLS_l1_penalty_one = JLS_l1_penalty_one,
                         chrs = chrs,
                         num_chr_parallel = 4)

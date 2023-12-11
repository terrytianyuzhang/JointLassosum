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
source("generate_synthetic_data_function.R")

JLS_population_weight_one <- 0.5 ###THIS NEED TO BE ONE OF THE CANDIDATE "GAMMA"
JLS_l1_penalty_one <- 0.005###THIS NEED TO BE ONE OF THE CANDIDATE "LAMBDA"

##the result of one JLS fit is the input of this pipeline
JLS_result_folder <- '/raid6/Tianyu/PRS/sharable/result/' ###POINT THE CODE TO THE OUTPUT OF JLS FITTING RESULT. AFTER RUNNING fit_JLS.R
large_population_type <- 'CEU'
small_population_type <- 'YRI'
synthetic_large_population_prefix <- '/raid6/Tianyu/PRS/sharable/data/CEU-chr21n22' ###THE SYNTHETIC POPULATION GENOTYPE DATA, ALL CHROMOSOMES
synthetic_small_population_prefix <- '/raid6/Tianyu/PRS/sharable/data/YRI-chr21n22'

synthetic_large_population_prefix_by_chr <- '/raid6/Tianyu/PRS/sharable/data/CEU-chr' ###THE SYNTHETIC POPULATION GENOTYPE DATA, SPLIT BY CHROMOSOMES
synthetic_small_population_prefix_by_chr <- '/raid6/Tianyu/PRS/sharable/data/YRI-chr'

large_population_GWAS_file <- '/raid6/Tianyu/PRS/sharable/data/large_population_GWAS_two_chr'#GWAS RESULTS, COPY THIS FROM fit_JLS.R
small_population_GWAS_file <- '/raid6/Tianyu/PRS/sharable/data/small_population_GWAS_two_chr'
large_population_GWAS_case_proportion <- 0.5
small_population_GWAS_case_proportion <- 0.5
chrs <- 21:22

###OUTPUT FOLDER
para_tuning_result_folder <- '/raid6/Tianyu/PRS/sharable_synthetic_tuning/result/' ###PLACE TO STORE PARAMETER TUNING RESULTS

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
                          GWAS_file = large_population_GWAS_file,
                          population_type = large_population_type,
                          para_tuning_result_folder = para_tuning_result_folder,
                          case_proportion = large_population_GWAS_case_proportion,
                          JLS_population_weight_one = JLS_population_weight_one,
                          JLS_l1_penalty_one = JLS_l1_penalty_one)
synthetic_label_given_PGS(JLS_result_folder = JLS_result_folder,
                          GWAS_file = small_population_GWAS_file,
                          population_type = small_population_type,
                          para_tuning_result_folder = para_tuning_result_folder,
                          case_proportion = small_population_GWAS_case_proportion,
                          JLS_population_weight_one = JLS_population_weight_one,
                          JLS_l1_penalty_one = JLS_l1_penalty_one)



#####STEP 3: SPLIT THE TRAINING AND VALIDATION SETS WITHIN THE SYNTHETIC POPULATION
set.seed(2019)
population_type <- 'CEU'
GWAS_file <- large_population_GWAS_file
synthetic_population_prefix_by_chr<- synthetic_large_population_prefix_by_chr

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

######SECTION 4: CALCULATE GWAS WITH SYNTHETIC DATA######


calculate_synthetic_GWAS(population_type = large_population_type,
                         para_tuning_result_folder = para_tuning_result_folder,
                         synthetic_population_prefix_by_chr = synthetic_large_population_prefix_by_chr,
                         GWAS_file = large_population_GWAS_file,
                         JLS_population_weight_one = JLS_population_weight_one,
                         JLS_l1_penalty_one = JLS_l1_penalty_one,
                         chrs = chrs,
                         num_chr_parallel = 1)
calculate_synthetic_GWAS(population_type = small_population_type,
                         para_tuning_result_folder = para_tuning_result_folder,
                         synthetic_population_prefix_by_chr = synthetic_small_population_prefix_by_chr,
                         GWAS_file = small_population_GWAS_file,
                         JLS_population_weight_one = JLS_population_weight_one,
                         JLS_l1_penalty_one = JLS_l1_penalty_one,
                         chrs = chrs,
                         num_chr_parallel = 8)











####

dir.create(paste0(ParameterTuningDirectory, '/Assoc/'),
           showWarnings = F,recursive = T)

####

for(population_index in 1:2){
  ####
  ancestry <- ancs[population_index]
  
  if(ancestry == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'
  }
  ####
  
  TrainSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", 
                                  ancestry, "-synthetic-train-index.Rdata")
  train.index <- get(load(TrainSampleIndexRFile))
  
  ####
  
  SyntheticYFile <- paste0(ParameterTuningDirectory, '/',
                           ancestry,'-SyntheticY', '.RData')
  SyntheticY <- get(load(SyntheticYFile))
  SyntheticY <- SyntheticY[train.index]
  
  ####
  
  for(chr in chrs){
    
    unlabeled_fam_file <- paste0(ParameterTuningDirectory,
                                 "/CHR/", ancestry, 
                                 "-chr", chr, "synthetic-train.fam")
    unlabeled_fam <- fread(unlabeled_fam_file)
    
    unlabeled_fam$V6 <- SyntheticY + 1
    fwrite(unlabeled_fam, file = unlabeled_fam_file, sep = ' ', col.names = FALSE)
    
    
    ####
    
    
    plink2.command=paste(plink2,"--nonfounders","--allow-no-sex",
                         "--bfile", paste0(ParameterTuningDirectory, "/CHR/", 
                                           ancestry,"-chr", chr,"synthetic-train"),
                         "--glm","allow-no-covars","omit-ref",
                         "--out", paste0(ParameterTuningDirectory, "/Assoc/", 
                                         ancestry,"-chr", chr,"synthetic-train"),
                         sep=" ")
    system(plink2.command)
    
  }
  ####
  
  all_GWAS <- data.table()
  for(chr in chrs){
    one_chr_GWAS <- fread(paste0(ParameterTuningDirectory, 
                               "/Assoc/", ancestry,"-chr", chr,
                               "synthetic-train.PHENO1.glm.logistic.hybrid"),
                        header = T)
    all_GWAS <- rbind(all_GWAS, one_chr_GWAS)
  }
  
  ####
  fwrite(all_GWAS, paste0(ParameterTuningDirectory, 
                          "/Assoc/", ancestry,"-",
                          "synthetic-train.PHENO1.glm.logistic.hybrid"),
         sep = ' ')
  
}

timing_file <- paste0(TimingResultDirectory,
                      '/timing_result.rds')
saveRDS(timing_result, file = timing_file)

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
synthetic_large_population_prefix <- '/raid6/Tianyu/PRS/sharable/data/CEU-chr21n22' ###THE SYNTHETIC POPULATION GENOTYPE DATA
synthetic_small_population_prefix <- '/raid6/Tianyu/PRS/sharable/data/YRI-chr21n22'
synthetic_large_population_prefix_by_chr <- '/raid6/Tianyu/PRS/sharable/data/CEU-chr' ###THE SYNTHETIC POPULATION GENOTYPE DATA
synthetic_small_population_prefix_by_chr <- '/raid6/Tianyu/PRS/sharable/data/YRI-chr'
large_population_GWAS_file <- '/raid6/Tianyu/PRS/sharable/data/large_population_GWAS_two_chr'#GWAS RESULTS, COPY THIS FROM fit_JLS.R
small_population_GWAS_file <- '/raid6/Tianyu/PRS/sharable/data/small_population_GWAS_two_chr'
large_population_GWAS_case_proportion <- 0.5
small_population_GWAS_case_proportion <- 0.5

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

#####STEP 3: SPLIT THE TRAINING AND VALIDATION SETS WITHIN THE SYNTHETIC POPULATION
set.seed(2019)
num_fold <- 5
chr <- 21
population_type <- 'YRI'
GWAS_file <- small_population_GWAS_file
synthetic_population_prefix_by_chr<- synthetic_small_population_prefix_by_chr
print('-----generate train and validation sample indices-----')

train_sample_index_r_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                     paste0(population_type, '_train_index_'),
                                                                     JLS_population_weight_one,
                                                                     JLS_l1_penalty_one,
                                                                     '.Rdata')
validation_sample_index_r_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                       paste0(population_type, '_validation_index_'),
                                                                       JLS_population_weight_one,
                                                                       JLS_l1_penalty_one,
                                                                       '.Rdata')
GWAS <- fread(GWAS_file, nrows = 3) ##I ONLY NEED THE SAMPLE SIZE INFORMATION
population_size <- as.numeric(GWAS$OBS_CT[1])

validation_index <- sort(sample(1:population_size, floor(population_size/num_fold)))
train_index <- (1:population_size)[-validation_index]

save(train_index, file = train_sample_index_r_file)
save(validation_index, file = validation_sample_index_r_file)
print('-----finished-----')

####
print('-----split the plink file-----')
train_sample_index_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                     paste0(population_type, '_train_index_'),
                                                                     JLS_population_weight_one,
                                                                     JLS_l1_penalty_one,
                                                                     '.txt')
validation_sample_index_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                          paste0(population_type, '_validation_index_'),
                                                                          JLS_population_weight_one,
                                                                          JLS_l1_penalty_one,
                                                                          '.txt')
synthetic_population_fam <- fread(paste0(synthetic_population_prefix, ".fam"))
train_fam <- synthetic_population_fam[train_index, c(1,2)] #only keep family id and withtin family id
fwrite(train_fam, train_sample_index_file, col.names = F, sep = " ")
val_fam <- synthetic_population_fam[validation_index, c(1,2)]
fwrite(val_fam, validation_sample_index_file, col.names = F, sep = " ")

train_genotype_folder_by_chr <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                          paste0(population_type, '_train_CHR_'),
                                                                          JLS_population_weight_one,
                                                                          JLS_l1_penalty_one,
                                                                          '/')
dir.create(train_genotype_folder_by_chr, showWarnings = F,recursive = T)

train_genotype_prefix <- paste0(train_genotype_folder_by_chr, "chr_", chr)

plink.command <- paste("plink --bfile", paste0(synthetic_population_prefix_by_chr, chr),
                       "--allow-no-sex",
                       "--keep", train_sample_index_file,
                       "--make-bed", "--out", train_genotype_prefix,
                       "--noweb", "--keep-allele-order",
                       sep = " ")
system(plink.command)

validation_genotype_folder_by_chr <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                          paste0(population_type, '_validation_CHR_'),
                                                                          JLS_population_weight_one,
                                                                          JLS_l1_penalty_one,
                                                                          '/')
dir.create(validation_genotype_folder_by_chr, showWarnings = F,recursive = T)

validation_genotype_prefix <- paste0(validation_genotype_folder_by_chr, "chr_", chr)

plink.command <- paste("plink --bfile", paste0(synthetic_population_prefix_by_chr, chr),
                       "--allow-no-sex",
                       "--keep", validation_sample_index_file,
                       "--make-bed", "--out", validation_genotype_prefix,
                       "--noweb", "--keep-allele-order",
                       sep = " ")
system(plink.command)

print(paste0('finished sample splitting for chr ', chr))
print('-----done-----')



mclapply(chrs, splitTrainValidation, anc = anc,
         train.index = train.index, val.index = val.index,
         ParameterTuningDirectory = ParameterTuningDirectory,
         TrainSampleIndexFile = TrainSampleIndexFile,
         ValidationSampleIndexFile = ValidationSampleIndexFile,
         plink = plink, mc.cores = 16, mc.preschedule = FALSE)
######SECTION 3: split training and validation individuals####
###### 1/(nfold) left out for validation ####
##########split training and validation########


for(i.set in 1:2){
  anc <- ancs[i.set]
  s.size <- s.sizes[i.set]
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'
  }
  
  TrainSampleIndexFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.txt")
  ValidationSampleIndexFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-validate-index.txt")
  TrainSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-train-index.Rdata")
  ValidationSampleIndexRFile <- paste0(ParameterTuningDirectory, "/", anc, "-synthetic-validate-index.Rdata")
  
  val.index <- sort(sample(1:s.size, floor(s.size/TrainTestNFold)))
  train.index <- (1:s.size)[-val.index] #this is in order
  
  save(val.index, file = ValidationSampleIndexRFile)
  save(train.index, file = TrainSampleIndexRFile)
  
  ####
  referece_panel_allSNP_name <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/", file.title)
  all_sample_psam <- fread(paste0(referece_panel_allSNP_name, ".psam"))
  train_psam <- all_sample_psam[train.index, c(1,2)] #only keep family id and withtin family id
  fwrite(train_psam, TrainSampleIndexFile, col.names = F, sep = " ")
  val_psam <- all_sample_psam[val.index, c(1,2)]
  fwrite(val_psam, ValidationSampleIndexFile, col.names = F, sep = " ")
  
  ####generate training and testing .fam files
  mclapply(chrs, splitTrainValidation, anc = anc,
          train.index = train.index, val.index = val.index,
          ParameterTuningDirectory = ParameterTuningDirectory,
          TrainSampleIndexFile = TrainSampleIndexFile,
          ValidationSampleIndexFile = ValidationSampleIndexFile,
          plink = plink, mc.cores = 16, mc.preschedule = FALSE)
  print(paste(anc, 'is done splitting'))
}

print('finished sample splitting')

######SECTION 4: calculate GWAS on synthetic data######
chrs <- 1:22
ancs <- c('CEU', 'YRI')

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

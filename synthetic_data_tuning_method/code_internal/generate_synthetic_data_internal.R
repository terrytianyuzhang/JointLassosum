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
JLS_l1_penalty_one <- 0.02236068###THIS NEED TO BE ONE OF THE CANDIDATE "LAMBDA"

##the result of one JLS fit is the input of this pipeline
JLS_result_folder <- '/raid6/Tianyu/PRS/sharable/result_internal/' ###POINT THE CODE TO THE OUTPUT OF JLS FITTING RESULT. AFTER RUNNING fit_JLS.R
large_population_type <- 'CEU'
small_population_type <- 'YRI'
large_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/CEU/CHR/CEU-chr'
small_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/YRI/CHR/YRI-chr'
large_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid'
small_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/YRI.TRN.PHENO1.glm.logistic.hybrid'
large_population_GWAS_case_proportion <- 0.5
small_population_GWAS_case_proportion <- 0.5



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
#
timing_result <- system.time({

  #####STEP 2: BASED ON THE PGS, ASSIGN EACH OF THEM A SYNTHETIC OUTCOME LABEL
  
  
risk.score.list <- get(load(paste0(ParameterTuningDirectory,
                                   '/riskscore.Rdata')))
#####we need to read figure out what is the original data noise level
ancs <- c('CEU', 'YRI')
CEUSampleSize <- sample_sizes$CEU$n.case + sample_sizes$CEU$n.control
YRISampleSize <- sample_sizes$YRI$n.case + sample_sizes$YRI$n.control

s.sizes <- c(CEUSampleSize, YRISampleSize)
caseProportion <- sample_sizes$CEU$n.case / CEUSampleSize ###the case proportion is the same for both populations

for(i.set in 1:2){
  TrainGWASFile <- paste0(work.dir, 'TRN/',
                          ancs[i.set],'.TRN.PHENO1.glm.logistic.hybrid')
  SyntheticYFile <- paste0(ParameterTuningDirectory, '/',
                           ancs[i.set],'-SyntheticY', '.RData')

  
}

print('generated synthetic outcome Y')


######SECTION 3: split training and validation individuals####
###### 1/(nfold) left out for validation ####
##########split training and validation########
set.seed(2019)

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
})
timing_file <- paste0(TimingResultDirectory,
                      '/timing_result.rds')
saveRDS(timing_result, file = timing_file)

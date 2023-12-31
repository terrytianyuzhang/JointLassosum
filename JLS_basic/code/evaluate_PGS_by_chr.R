gc()

#####THIS PIECE OF CODE USES RESULTS GENERATED BY fit_JLS.R

library(data.table)
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(R.utils)
# library(pROC) # for AUC 
# library(pryr) # check memory useage
JLS_population_weight <- c(0, 0.5, 1) #gamma parameter in the paper
JLS_population_weight_one <- JLS_population_weight[2]

JLS_result_prefix <- '/raid6/Tianyu/PRS/sharable/result/JLS_result_weight_is'

###THIS PIECE OF CODE ALLOW THE PLINK FILES ARE STORED BY CHROMOSOMES
small_population_reference_prefix_by_chr<- '/raid6/Tianyu/PRS/sharable/data/YRI-chr'
chrs <- 21:22

###THE JLS MODELS ARE STORED SEPARATELY FOR DIFFERENT POPULATION WEIGHTS
JLS_regression_coefficient_file <- paste0(JLS_result_prefix, 
                                          sprintf("%.2f",JLS_population_weight_one), '_coefficient.txt')
JLS_regression_coefficient <- fread(JLS_regression_coefficient_file)

JLS_result_one_weight_file <- paste0(JLS_result_prefix, 
                                     sprintf("%.2f",JLS_population_weight_one), '.Rdata')

JLS_result_one_weight <- get(load(JLS_result_one_weight_file))

####NOW WE HAVE ALL THE COEFFICIENTS (BETA) FOR ONE WEIGHT (GAMMA)
####EACH OF THE COEFFICIENT VECTOR CORRESPONDS TO ONE L1 PENALTY (LAMBDA)
for(l1_penalty_index in 1:length(JLS_result_one_weight$lambda)){

  ####WE NEED TO STORE THE BETA IN CERTAIN FORMAT TO LEVERAGE THE FAST PLINK SOFTWARE FOR PGS EVALUATION
  effect_size_df <- data.frame(SNP = JLS_regression_coefficient$ID, 
                               A1 = unlist(lapply(strsplit(JLS_regression_coefficient$ID, ":"),`[[`,4)),
                               BETA = JLS_regression_coefficient[, ..l1_penalty_index])
  
  selected_coefficient_file <- paste0(JLS_result_prefix, 
                                      sprintf("%.2f",JLS_population_weight_one),
                                      '_l1_penalty_is_',
                                      sprintf("%.4f",JLS_result_one_weight$lambda[l1_penalty_index]), '_coefficient.txt')
  write.table(effect_size_df, selected_coefficient_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  evaluate_PGS_by_chr_parallel(JLS_result_prefix = JLS_result_prefix,
                               JLS_population_weight = JLS_population_weight_one,
                               JLS_l1_penalty = JLS_result_one_weight$lambda[l1_penalty_index],
                               genotype_plink_file = small_population_reference_prefix_by_chr,
                               chrs = chrs)
}

evaluate_PGS_by_chr_parallel <- function(JLS_result_prefix,
                                         JLS_population_weight,
                                         JLS_l1_penalty,
                                         genotype_plink_file, 
                                         chrs,
                                         num_chr_para = 11){
  calculate_PGS_one_chr <- function(chr){
    PGS_file <- paste0(JLS_result_prefix, 
                       sprintf("%.2f",JLS_population_weight),
                       '_l1_penalty_is_',
                       sprintf("%.4f",JLS_l1_penalty), '_PGS_chr_',chr)
    
    ###NOW THE REAL WORK IS HAPPENING
    plink2.command = paste("plink2 --nonfounders","--allow-no-sex","--threads", 2,"--memory", 25000,
                           "--bfile", paste0(genotype_plink_file, chr),
                           "--score", selected_coefficient_file, "header-read",1,2,
                           "--score-col-nums",3,
                           "--out", PGS_file,
                           sep=" ")
    
    system(plink2.command)
  }  
  
  mclapply(chrs,
           FUN = calculate_PGS_one_chr,
           mc.cores = num_chr_para, 
           mc.preschedule = F, 
           mc.silent = F)
  
  NAMED_ALLELE_DOSAGE_SUM <- ALLELE_CT <- 0 #TOTAL ALLELE
  for(chr in chrs){
    PGS_file <- paste0(JLS_result_prefix, 
                       sprintf("%.2f",JLS_population_weight),
                       '_l1_penalty_is_',
                       sprintf("%.4f",JLS_l1_penalty), '_PGS_chr_',chr, '.sscore')
    PGS <- fread(PGS_file, header = T)
    ALLELE_CT_chr <- as.numeric(PGS[1,]$ALLELE_CT)
    ALLELE_CT <- ALLELE_CT + ALLELE_CT_chr
    
    if(chr == chrs[1]){
      individual_PGS <- ALLELE_CT_chr * PGS[,6]
      individual_NAMED_ALLELE_DOSAGE_SUM <- PGS[,5]
    }else{
      individual_PGS <- individual_PGS + ALLELE_CT_chr * PGS[,6]
      individual_NAMED_ALLELE_DOSAGE_SUM <- individual_NAMED_ALLELE_DOSAGE_SUM + PGS[,5]
    }
  }
  individual_PGS <- individual_PGS/ALLELE_CT
  
  PGS_overall <- PGS
  PGS_overall$ALLELE_CT <- ALLELE_CT
  PGS_overall$NAMED_ALLELE_DOSAGE_SUM <- individual_NAMED_ALLELE_DOSAGE_SUM
  PGS_overall[,6] <- individual_PGS
  
  PGS_file <- paste0(JLS_result_prefix, 
                     sprintf("%.2f",JLS_population_weight),
                     '_l1_penalty_is_',
                     sprintf("%.4f",JLS_l1_penalty), '_PGS_all_chr', '.sscore')
  fwrite(PGS_overall, file = PGS_file, sep = "\t")
}



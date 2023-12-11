
file_name_generator_weight_and_l1_penalty <- function(folder, prefix, weight, l1_penalty = NA, suffix){
  first_half <- paste0(folder, prefix, "weight_is_", sprintf("%.2f",weight))
  if(is.na(l1_penalty)){
    file_name <- paste0(first_half, suffix)
  }else{
    file_name <- paste0(first_half, '_l1_penalty_is_', sprintf("%.4f",l1_penalty), suffix)
  }
  return(file_name)
}

####CALCULATE THE RISK SCORE OF SOME SYNTHETIC INDIVIDUALS
####THE GENOTYPE INFORMATION DOES NOT NEED TO CORRESPONDS TO REAL INDIVIDUALS
####THE REGRESSION COEFFICIENT VECTOR IS CALCULATED FROM THE PRELIMINARY FIT OF JLS

predict_PGS_given_coefficient <- function(JLS_result_folder, para_tuning_result_folder,
                                          population_type, synthetic_population_prefix,
                                          JLS_population_weight_one, JLS_l1_penalty_one){
  
  print("-----now loading the JLS model-----")
  JLS_regression_coefficient_file <- paste0(JLS_result_folder, 
                                           'JLS_result_weight_is',
                                           sprintf("%.2f",JLS_population_weight_one),
                                           '_coefficient.txt')
  JLS_regression_coefficient <- fread(JLS_regression_coefficient_file)
  
  
  JLS_result_one_weight_file <- paste0(JLS_result_folder, 
                                       'JLS_result_weight_is',
                                       sprintf("%.2f",JLS_population_weight_one),
                                       '.Rdata')
  JLS_result_one_weight <- get(load(JLS_result_one_weight_file))
  
  print("-----finished-----")
  
  print("-----save the regression coefficient in the required format by plink-----")
  l1_penalty_index <- which(abs(JLS_result_one_weight$lambda - JLS_l1_penalty_one) < 10^(-5))[1]
  effect_size_df <- data.frame(SNP = JLS_regression_coefficient$ID, 
                               A1 = unlist(lapply(strsplit(JLS_regression_coefficient$ID, ":"),`[[`,4)),
                               BETA = JLS_regression_coefficient[, ..l1_penalty_index])
  
  selected_coefficient_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                         'JLS_result_',
                                                                         JLS_population_weight_one,
                                                                         JLS_l1_penalty_one,
                                                                         '_coefficient.txt')
  write.table(effect_size_df, selected_coefficient_file, sep = "\t", quote = FALSE, row.names = FALSE)
  print("-----done-----")
  
  ####PLACE TO STORE THE PGS RESULTS
  PGS_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                        paste0(population_type, '_prediction_score_'),
                                                        JLS_population_weight_one,
                                                        JLS_l1_penalty_one,
                                                        NULL)
  
  print("-----push the PGS prediction work to plink-----")
  ###NOW THE REAL WORK IS HAPPENING
  plink2.command <- paste("plink2 --nonfounders","--allow-no-sex","--threads", 8,"--memory", 25000,
                          "--bfile", synthetic_population_prefix ,
                          "--score", selected_coefficient_file, "header-read",1,2,
                          "--score-col-nums",3,
                          "--out",PGS_file,
                          sep=" ")
  
  system(plink2.command)
  print("-----plink finished-----")
}

####GENERATE SYNTHETIC LABELS
####SAMPLE SIZE, CASE-CONTROL RATIO SHOULD BE SIMILAR TO THAT OF THE ORIGINAL GWAS STUDY

#UNCOMMENT FOR TESTING PURPOSE
# GWAS_file = small_population_GWAS_file
# population_type = small_population_type
# case_proportion = small_population_GWAS_case_proportion
# model_based_formula = FALSE
synthetic_label_given_PGS <- function(JLS_result_folder, GWAS_file,
                                      population_type, para_tuning_result_folder,
                                      case_proportion,   ###case_proportion > 0.5 means there is more case than control 
                                      JLS_population_weight_one,
                                      JLS_l1_penalty_one, extra_scaling = 1,
                                      model_based_formula = FALSE){
  
  print('-----load the risk score information generated last step-----')
  PGS_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                        paste0(population_type, '_prediction_score_'),
                                                        JLS_population_weight_one,
                                                        JLS_l1_penalty_one,
                                                        '.sscore')
  ###READ IN THE PREDICTED PGS SCORE FOR EACH INDIVIDUAL IN THE SYNTHETIC POPULATION
  PGS_complete <- fread(PGS_file, header = T)
  risk_score <- as.matrix(PGS_complete[, 6])
  risk_score <- risk_score - mean(risk_score)
  
  print('-----done-----')
  
  ##load the original GWAS results
  GWAS <- fread(GWAS_file, nrows = 3) ##I ONLY NEED THE SAMPLE SIZE INFORMATION
  population_size <- as.numeric(GWAS$OBS_CT[1])
  
  if(model_based_formula == TRUE){
    ###model-based formula, derivation included in the paper
    
    GWAS <- fread(GWAS_file) 
    
    print('-----load the preliminary regression coefficients-----')
    ###READ IN THE PRELIMINARY BETA
    JLS_result_one_weight_file <- paste0(JLS_result_folder, 
                                         'JLS_result_weight_is',
                                         sprintf("%.2f",JLS_population_weight_one),
                                         '.Rdata')
    JLS_result_one_weight <- get(load(JLS_result_one_weight_file))
    
    l1_penalty_index <- which(abs(JLS_result_one_weight$lambda - JLS_l1_penalty_one) < 10^(-15))
    beta_preliminary <- as.matrix(JLS_result_one_weight$beta[, l1_penalty_index])
    print('-----done-----')
    
    pheno_gene_correlation <- p2cor(p = GWAS$P, n = population_size, sign = log(GWAS$OR)) #original correlation vector
    design_factor <- sqrt(case_proportion * (1 - case_proportion)) 
    rescaled_variance <- sum(risk_score^2) #this is the (rescaled) variance of the risk scores
    slope <- pheno_gene_correlation %*% beta_preliminary * sqrt(population_size) * design_factor / rescaled_variance 
    slope <- as.numeric(slope) * extra_scaling ###sometimes may need a extra.scaling != 1 to generate better cv data sets
    
    case_probability <-  slope * risk_score + case_proportion
  }else{
    risk_score <- risk_score/ max(max(risk_score), -min(risk_score))
    case_probability <-  min(case_proportion, 1-case_proportion) * risk_score + case_proportion
  }
  
  
  #make sure mu is between 0 and 1 since it is a probability
  # print(paste0('number of beyond 0/1 range', sum(case_probability > 1 | case_probability<0)))
  case_probability[case_probability > 1] <- 1
  case_probability[case_probability < 0] <- 0
  
  ###DRAW THE SYNTEHTIC LABEL
  synthetic_label <- rbinom(population_size, 1, case_probability)
  print('-----summary of the synthetic label-----')
  print(summary(synthetic_label))
  
  synthetic_label_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                    paste0(population_type, '_synthetic_label_'),
                                                                    JLS_population_weight_one,
                                                                    JLS_l1_penalty_one,
                                                                    '.Rdata')
  save(synthetic_label, file = synthetic_label_file)
}

split_train_validation <- function(chr, population_type,
                                   para_tuning_result_folder,
                                   synthetic_population_prefix_by_chr,
                                   GWAS_file, 
                                   JLS_population_weight_one,
                                   JLS_l1_penalty_one,
                                   num_fold = 5){
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
  synthetic_population_fam <- fread(paste0(synthetic_population_prefix_by_chr, chr, ".fam"))
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
}

calculate_synthetic_GWAS <- function(population_type,
                                     para_tuning_result_folder,
                                     synthetic_population_prefix_by_chr,
                                     GWAS_file, 
                                     JLS_population_weight_one,
                                     JLS_l1_penalty_one,
                                     chrs,
                                     num_chr_parallel = 1){
  train_sample_index_r_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                         paste0(population_type, '_train_index_'),
                                                                         JLS_population_weight_one,
                                                                         JLS_l1_penalty_one,
                                                                         '.Rdata')
  train_index <- get(load(train_sample_index_r_file))
  
  ####TRAINING LABEL
  synthetic_label_file <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                    paste0(population_type, '_synthetic_label_'),
                                                                    JLS_population_weight_one,
                                                                    JLS_l1_penalty_one,
                                                                    '.Rdata')
  synthetic_label <- get(load(synthetic_label_file))
  synthetic_label <- synthetic_label[train_index]
  
  GWAS_folder <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                           paste0(population_type, '_synthetic_GWAS_'),
                                                           JLS_population_weight_one,
                                                           JLS_l1_penalty_one,
                                                           '/')
  
  calculate_GWAS_one_chr <- function(chr){
    ####MODIFY THE PHENO COLUMN TO BE 1&2 INSTEAD OF 0&1, NOT SURE IF NECESSARY
    train_genotype_folder_by_chr <- file_name_generator_weight_and_l1_penalty(para_tuning_result_folder, 
                                                                              paste0(population_type, '_train_CHR_'),
                                                                              JLS_population_weight_one,
                                                                              JLS_l1_penalty_one,
                                                                              '/')
    
    train_genotype_prefix <- paste0(train_genotype_folder_by_chr, "chr_", chr)
    
    unlabeled_fam_file <- paste0(train_genotype_prefix, '.fam')
    unlabeled_fam <- fread(unlabeled_fam_file)
    
    unlabeled_fam$V6 <- synthetic_label + 1
    fwrite(unlabeled_fam, file = unlabeled_fam_file, sep = ' ', col.names = FALSE)
    
    
    dir.create(GWAS_folder, showWarnings = F,recursive = T)
    GWAS_file <- paste0(GWAS_folder, 'gwas_chr_', chr)
    
    plink2_command <- paste("plink2 --nonfounders","--allow-no-sex",
                            "--bfile", train_genotype_prefix,
                            "--glm","allow-no-covars","omit-ref",
                            "--out", GWAS_file,
                            sep=" ")
    system(plink2_command)
  }
  
  mclapply(chrs, calculate_GWAS_one_chr, mc.cores = num_chr_parallel, mc.preschedule = FALSE)
  
  all_GWAS <- data.table()
  for(chr in chrs){
    GWAS_file <- paste0(GWAS_folder, 'gwas_chr_', chr, '.PHENO1.glm.logistic.hybrid')
    one_chr_GWAS <- fread(GWAS_file, header = T)
    all_GWAS <- rbind(all_GWAS, one_chr_GWAS)
  }
  
  ####
  all_GWAS_file <- paste0(GWAS_folder, 'gwas_all.PHENO1.glm.logistic.hybrid')
  fwrite(all_GWAS, all_GWAS_file, sep = ' ')
}

splitTrainValidation <- function(chr = NULL, anc, train.index, val.index, 
                            ParameterTuningDirectory, 
                            TrainSampleIndexFile,
                            ValidationSampleIndexFile,
                            plink){
  
  if(!is.null(chr)){
    print(paste0('splitting train and validation for chr ', chr))
  }else{
    print(paste0('splitting train and validation for all the SNPs'))
  }
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  if(!is.null(chr)){
    ##this is the reference panel genotype data befpre splitting
    referece.panel.name <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/", file.title,"/CHR/",file.title,"-chr", chr)
    # boost.data.folder <- paste0("/raid6/Tianyu/PRS/BootData/", file.title,"/CHR")
    # boost.data.train.index <- paste0(boost.data.folder, "/boost-train-index.txt")
    # boost.data.val.index <- paste0(boost.data.folder, "/boost-val-index.txt")
    # TrainSampleDataPrefix <- paste0(boost.data.folder, "/", file.title,"-chr", chr, "synthetic-train")
    # ValidateSampleDataPrefix <- paste0(boost.data.folder, "/", file.title,"-chr", chr, "synthetic-val")
    # 
    dir.create(paste0(ParameterTuningDirectory, "/CHR/"),
               showWarnings = F,recursive = T)
    TrainSampleDataPrefix <- paste0(ParameterTuningDirectory, "/CHR/", 
                                    anc, "-chr", chr, "synthetic-train")
    ValidateSampleDataPrefix <- paste0(ParameterTuningDirectory, "/CHR/", 
                                       anc, "-chr", chr, "synthetic-val")
    
    #complete panel information
    # full.fam <- fread(paste0("/raid6/Tianyu/PRS/bert_sample/",anc,".TUNE/CHR/",anc,".TUNE-chr",chr,".fam"))
    # full.fam <- fread(paste0(referece.panel.name, ".fam"))
    
    # train.fam <- full.fam[train.index, c(1,2)] #only keep family id and withtin family id
    # #write down which are in the training set
    # # fwrite(train.fam, paste0("/raid6/Tianyu/PRS/BootData/",anc,".TUNE/CHR/",anc,".TUNE-boost-train-index.txt"),
    # #        col.names = F, sep = " ")
    # fwrite(train.fam, TrainSampleIndexFile,
    #        col.names = F, sep = " ")
    
    #split the reference panel into training and validation sets
    plink.command <- paste(plink, "--bfile", referece.panel.name,
                           "--allow-no-sex",
                           "--keep", TrainSampleIndexFile,
                           "--make-bed", "--out", TrainSampleDataPrefix,
                           "--noweb", "--keep-allele-order",
                           sep = " ")
    system(plink.command)
    
    # val.fam <- full.fam[val.index, c(1,2)]
    # #write down which are in the validation set
    # fwrite(val.fam, ValidationSampleIndexFile,
    #        col.names = F, sep = " ")
    
    plink.command <- paste(plink, "--bfile", referece.panel.name,
                           "--allow-no-sex",
                           "--keep", ValidationSampleIndexFile,
                           "--make-bed", "--out", ValidateSampleDataPrefix,
                           "--noweb", "--keep-allele-order",
                           sep = " ")
    system(plink.command)
    
    print(paste0('finished sample splitting for chr ', chr))
  }else{
    ##this is the reference panel genotype data befpre splitting
    referece.panel.name <- paste0("/raid6/Tianyu/PRS/bert_sample/ReferencePopulation-Package/",file.title, "/", file.title)
    
    TrainSampleDataPrefix <- paste0(ParameterTuningDirectory, "/", 
                                    anc, "synthetic-train")
    ValidateSampleDataPrefix <- paste0(ParameterTuningDirectory, "/", 
                                       anc, "synthetic-val")
    
    #complete panel information
    full.fam <- fread(paste0(referece.panel.name, ".psam"))
    
    train.fam <- full.fam[train.index, c(1,2)] #only keep family id and withtin family id
    fwrite(train.fam, TrainSampleIndexFile,
           col.names = F, sep = " ")
    
    #split the reference panel into training and validation sets
    plink.command <- paste(plink, "--bfile", referece.panel.name,
                           "--allow-no-sex",
                           "--keep", TrainSampleIndexFile,
                           "--make-bed", "--out", TrainSampleDataPrefix,
                           "--noweb", "--keep-allele-order",
                           sep = " ")
    system(plink.command)
    
    val.fam <- full.fam[val.index, c(1,2)]
    #write down which are in the validation set
    fwrite(val.fam, ValidationSampleIndexFile,
           col.names = F, sep = " ")
    
    plink.command <- paste(plink, "--bfile", referece.panel.name,
                           "--allow-no-sex",
                           "--keep", ValidationSampleIndexFile,
                           "--make-bed", "--out", ValidateSampleDataPrefix,
                           "--noweb", "--keep-allele-order",
                           sep = " ")
    system(plink.command)
    
    print(paste0('finished sample splitting for all the SNPs'))
  }
}

PairwiseCorrelationSyntheticData <- function(chr, anc, SyntheticY, ParameterTuningDirectory){
  print(paste0('chr is ', chr))
  ##calculate pairwise correlation
  
  # chr_loc <- as.numeric(gsub(':.*$','',names(beta0)))
  # snp <- names(beta0)[chr_loc == chr]
  
  if(anc == 'CEU'){
    file.title <- 'CEU-20K'
  }else{
    file.title <- 'YRI-4K'  
  }
  
  # boost.data.folder <- paste0("/raid6/Tianyu/PRS/BootData/", file.title,"/CHR")
  # boost.data.train.data <- paste0(boost.data.folder, "/", file.title,"-chr", chr, "boost-train")
  TrainSampleDataPrefix <- paste0(ParameterTuningDirectory, "/",
                                  anc, "-chr", chr, "synthetic-train")
  ####load training genotype data
  system.time(gnt<-read.plink(bed=paste0(TrainSampleDataPrefix,".bed"),
                              bim=paste0(TrainSampleDataPrefix,".bim"),
                              fam=paste0(TrainSampleDataPrefix,".fam"))
  )
  # system.time(gnt<-read.plink(bed=paste0(boost.data.train.data,".bed"),
  #                             bim=paste0(boost.data.train.data,".bim"),
  #                             fam=paste0(boost.data.train.data,".fam"))
  # )
  
  # transform the genotypes to a matrix, and reverse the call count. This snpStats package counts the A2 allele, not the A1
  system.time(gnt<-2-as(gnt$genotypes,Class="numeric"))
  
  print(paste0("the size of chr ", chr, " genotype matrix is (GB)"))
  print(object.size(gnt)/1e9)
  # center the columns
  system.time(gnt <-gnt - rep(1, nrow(gnt)) %*% t(colMeans(gnt)))
  # normalize the calls to the unit 1 norm
  system.time(gnt<-normalize.cols(gnt,method="euclidean",p=2))
  
  ###calculate pairwise correlation
  print('calculating pairwise correlation')
  SNP_number <- NCOL(gnt)
  sample_size <- NROW(gnt)
  
  
  boot.cor <- matrix(0, nrow = SNP_number, ncol = 1)
  for(j in 1:SNP_number){
    ###after the proper normalization, correlation is the same as inner product
    boot.cor[j,1] <- crossprod(SyntheticY[,1], gnt[,j] )
  }
  rm(gnt)
  #####make some "fake" GWAS results
  print('translate correlation into p-values')
  # SyntheticGWASOneChromosome <- list()
  
  # for(i in 1:B){
  t_stat <- boot.cor[,1] * sqrt(sample_size - 2)/sqrt(1 - (boot.cor[,1])^2)
  SyntheticGWASOneChromosome <- data.frame(P = 2*pt(q = abs(t_stat), 
                                           df = sample_size - 2, lower.tail = F), #two-sided p-value
                                           OR = 2*as.numeric(boot.cor[,1]>0) + 0.5*as.numeric(boot.cor[,1]<0))
  SyntheticGWASOneChromosome$n <- sample_size
  #   SyntheticGWASList[[i]] <- fake_GWAS
  # }
  
  return(SyntheticGWASOneChromosome)
}

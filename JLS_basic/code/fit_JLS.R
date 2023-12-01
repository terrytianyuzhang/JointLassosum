gc()
options(stringsAsFactors = F)

if(!exists("i.sim")){
  i.sim <- 609
}

#### load the functions that are needed
# source("simulation-functions.R")
# source("general_pipeline_parameters.R")
source("Lassosum_function/parseselect.R")
source("Lassosum_function/parseblocks.R")
source("Lassosum_function/ncol.bfile.R")
source("Lassosum_function/nrow.bfile.R")
source("Lassosum_function/read.table2.R")
source("Lassosum_function/selectregion.R")
source("Lassosum_function/parse.pheno.covar.R")
source("Lassosum_function/myelnet.R")
source("Lassosum_function/mylassosum.R")
source("Lassosum_function/splitgenome.R")
source("Lassosum_function/validation.R")
source("Lassosum_function/merge.mylassosum.R")
Rcpp::sourceCpp("Lassosum_function/myfunctions.cpp")

library(data.table)
library(pryr) # check memory useage
library(lassosum) #transform p value to correlation
library(doParallel) # foreach
library(pROC) # for AUC 
library(R.utils)

source('JLS_function.R')

###INPUT FILE LOCATION
# large_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/CEU.TRN.PHENO1.glm.logistic.hybrid'
# small_population_GWAS_file <- '/raid6/Ron/prs/data/bert_sample/YRI.TRN.PHENO1.glm.logistic.hybrid'
# large_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/CEU/CHR/CEU-chr'
# small_population_reference_prefix <- '/raid6/Tianyu/PRS/SimulationPipeline/Data/Reference-LDblocks/YRI/CHR/YRI-chr'

large_population_GWAS_file <- '/raid6/Tianyu/PRS/sharable/data/large_population_GWAS_two_chr'
small_population_GWAS_file <- '/raid6/Tianyu/PRS/sharable/data/small_population_GWAS_two_chr'
large_population_reference_prefix <- '/raid6/Tianyu/PRS/sharable/data/CEU-chr'
small_population_reference_prefix <- '/raid6/Tianyu/PRS/sharable/data/YRI-chr'

# large_population_GWAS_file <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/large_population_GWAS_two_chr'
# small_population_GWAS_file <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/small_population_GWAS_two_chr'
# large_population_reference_prefix <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/CEU-chr'
# small_population_reference_prefix <- '/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/data/YRI-chr'
# gene_ID_chromosome_file <- '/raid6/Ron/prs/data/bert_sample/CEU.TRN/CEU.TRN.pvar'


###OTHER METADATA
large_population_type <- 'CEU'
small_population_type <- 'YRI'
# small_population_type <- 'ASN'

###HYPERPARAMETER CANDIDATES
JLS_population_weight <- c(0.5, 0.8)
# JLS_l1_penalty <- exp(seq(log(0.007), log(0.05), length.out=5))
JLS_l1_penalty <- exp(seq(log(0.1), log(0.05), length.out = 2))
JLS_shrinkage <- c(0.75)
chromosome <- 21:22
mem.limit <- 2e10

###GIVEN THE ABOVE INFORMATION, FIT THE MODEL FOR ONE TIME
JLS_population_weight_one <- JLS_population_weight[1]

####LOAD THE LD BLOCK BOUNDARY INFORMATION
map_population_type_to_Berisa_label <- function(population_type){
  if(population_type %in% c("CEU", "EUR")){
    Berisa_label <- "EUR"
  }else if(population_type %in% c("YRI", "AFR")){
    Berisa_label <- "AFR"
  }else{
    print(paste0('I cannot find a Berisa label for ', population_type))
  }
  return(Berisa_label)
}

LD_block_boundary <- list()
large_population_LD_file <- paste0("data/Berisa.", 
                                   map_population_type_to_Berisa_label(large_population_type), ".hg38.bed")
LD_block_boundary[[large_population_type]] <- read.table2(system.file(large_population_LD_file, 
                                                                      package="lassosum"), 
                                                          header=T)
small_population_LD_file <- paste0("data/Berisa.", 
                                   map_population_type_to_Berisa_label(small_population_type), ".hg38.bed")
LD_block_boundary[[small_population_type]] <- read.table2(system.file(small_population_LD_file, 
                                                                      package="lassosum"), 
                                                          header=T)

###NOW WE ARRANGE THE REFERENCE PANEL FILE NAMES FOR PARALLELIZATION (OVER THE CHROMOSOMES)
reference_file <- list()
for(chromosome_one in chromosome){
  reference_file[[chromosome_one]] <- list()
  reference_file[[chromosome_one]][[large_population_type]] <- paste0(large_population_reference_prefix, chromosome_one)
  reference_file[[chromosome_one]][[small_population_type]] <- paste0(small_population_reference_prefix, chromosome_one)
}


####LOAD THE GWAS RESULTS AND TRANSLATE THEM INTO CORRELATION INFORMATION
print('-----now loading the GWAS results-----')
# gene_ID_chromosome <- fread(gene_ID_chromosome_file, header=T, data.table=F)
large_population_GWAS <-  fread(large_population_GWAS_file, header = T, data.table = F)
large_population_size <- as.numeric(large_population_GWAS$OBS_CT[1])
pheno_gene_correlation <- data.frame(CHR = large_population_GWAS$`#CHROM`,
                                     ID = large_population_GWAS$ID)
pheno_gene_correlation[,large_population_type] <- p2cor(p = large_population_GWAS$P, 
                                                        n = large_population_size, 
                                                        sign=log(large_population_GWAS$OR))

small_population_GWAS <-  fread(small_population_GWAS_file, header = T, data.table = F)
small_population_size <- as.numeric(small_population_GWAS$OBS_CT[1])
pheno_gene_correlation[,small_population_type] <- p2cor(p = small_population_GWAS$P, 
                                                        n = small_population_size, 
                                                        sign=log(small_population_GWAS$OR))
rownames(pheno_gene_correlation) <- pheno_gene_correlation$ID
print('-----finished-----')

######NOW START MODEL FITTING
JLS_result_by_chr <- mclapply(chromosome,
                              FUN = split_chromosome_population_n_send_to_worker,
                              JLS_population_weight_one = JLS_population_weight_one,
                              JLS_l1_penalty = JLS_l1_penalty,
                              JLS_shrinkage = JLS_shrinkage,
                              large_population_type = large_population_type,
                              small_population_type = small_population_type,
                              LD_block_boundary = LD_block_boundary,
                              reference_file = reference_file,
                              pheno_gene_correlation = pheno_gene_correlation,
                              mem.limit = mem.limit,
                              mc.cores = mymc.cores,
                              mc.preschedule = F)
# large_population_GWAS_two_chr <- large_population_GWAS[large_population_GWAS$`#CHROM` %in% c(21,22), ]
# small_population_GWAS_two_chr <- small_population_GWAS[small_population_GWAS$`#CHROM` %in% c(21,22), ]
# large_population_GWAS_two_chr_file <- '/raid6/Tianyu/PRS/sharable/large_population_GWAS_two_chr'
# small_population_GWAS_two_chr_file <- '/raid6/Tianyu/PRS/sharable/small_population_GWAS_two_chr'
# fwrite(large_population_GWAS_two_chr, file = large_population_GWAS_two_chr_file)
# fwrite(small_population_GWAS_two_chr, file = small_population_GWAS_two_chr_file)


###PARALLEL




# setting.title <- 'CEU1aYRI2a22Chr'
##### general setup
# i.sim="NC8000"
# load(paste0("Work/Sim-C0.80-Y0.60-rep",i.sim,"/simulation-params.RData"))
# load("/raid6/Ron/prs/data/bert_sample/simulation-params.RData")

# set directories for this simulation
# main.dir=params$run.info$main.dir
# work.dir=params$run.info$work.dir

# load the parameters for this simulation
# load(paste0("Work/Sim-",i.sim,"/simulation-params.RData"))
# main.dir=params$run.info$main.dir #"/raid6/Tianyu/PRS/SimulationPipeline/"
# work.dir=params$run.info$work.dir #"/raid6/Tianyu/PRS/SimulationPipeline/Work/Sim-800/"

### ancestries
gwasANC=c("CEU","YRI")
# #### set the parameters
# GAMMA = c(0.2, 0.5, 0.8)
# ###!!!!!
# lambda=exp(seq(log(0.0025), log(0.025), length.out=10))
# lambda[1] <- 0.00025
# lambda=exp(seq(log(0.007), log(0.05), length.out=10))  ##try larger lambda

# lambda=exp(seq(log(0.001), log(0.025), length.out=10))
# shrink=.9

####sample sizes
# N1 <- params$CEU.TRN$n.case + params$CEU.TRN$n.control
# N2 <- params$YRI.TRN$n.case + params$YRI.TRN$n.control
N1 <- sample_sizes$CEU$n.case + sample_sizes$CEU$n.control
N2 <- sample_sizes$YRI$n.case + sample_sizes$YRI$n.control
input.df <- data.frame(gamma = GAMMA,
                       N1 = rep(N1, length(GAMMA)),
                       N2 = rep(N2, length(GAMMA)))
###!!!!this is not correct
# N=seq(20000,20000,4000)
# input.df=data.frame(N1=rep(N,each=length(N)),N2=rep(N,length(N)))
# input.df=data.frame(gamma=rep(GAMMA,each=nrow(input.df)),N1=rep(input.df$N1,length(GAMMA)),N2=rep(input.df$N2,length(GAMMA)))
# 
# ###this is the correct sample size
# input.df <- data.frame(gamma = GAMMA, N1 = 20000, N2 = 4000)
# memory limit
mem.limit <- 2*10e9

# end of the wrapper function

system.time(re.wrapper<-mclapply(1:nrow(input.df),wrapperFunction,input.df=input.df,
                                 gwasANC = gwasANC, lambda=lambda, shrink=shrink,
                                 main.dir=main.dir,work.dir=work.dir,CHR= 1:22, mem.limit=2e10,
                                 mc.cores=3, mc.preschedule = F, mc.silent=F))

rm.list=ls()
rm.list=rm.list[!(rm.list %in% c("i.sim","sims", "plink", "plink2"))]
rm(list = rm.list); flush.console()



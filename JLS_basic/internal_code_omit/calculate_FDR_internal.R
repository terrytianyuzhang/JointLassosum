library(data.table)
library(ggplot2)
library(ggpubr)


# plot_list <- list()
fairness_data <- data.frame()
setting_index <- 1

rawSimulationAUC <- readRDS(paste0('/Users/tianyuzhang/Documents/GitHub/JointLassosum/JLS_basic/internal_data_omit/YRI-PGS-overlap-100_2023-03-31.RDS'))
testing_population <- 'YRI'
repeat_index <- 1
repitition_name <- paste0('rep.', repeat_index)
raw_PRS_onesimulation <- data.table(rawSimulationAUC[[testing_population]][[repitition_name]])
raw_PRS_onesimulation <- melt(raw_PRS_onesimulation, 
                              id.vars = c("ID", "PHENO"), 
                              variable.name = "method", 
                              value.name = "pgs_score")
raw_PRS_onesimulation <- raw_PRS_onesimulation[method == 'JLS',]

print('----please formulate your data in the format of raw_PRS_onesimulation-----')
print('----PHENO == 1 corresponds to control individual-----')
print(raw_PRS_onesimulation)
# > print(raw_PRS_onesimulation)
# ID PHENO method  pgs_score
# 1:     YRI.78     2    JLS  1.6350823
# 2:    YRI.128     1    JLS -1.1602135
# 3:    YRI.132     2    JLS  0.3089623
# 4:    YRI.194     2    JLS  0.1308629
# 5:    YRI.226     2    JLS -0.5851154
# ---                                   
#   3996: YRI.199827     2    JLS  1.0477648
# 3997: YRI.199910     1    JLS -0.4124180
# 3998: YRI.199970     1    JLS -2.2842603
# 3999: YRI.199971     2    JLS  1.1195098
# 4000: YRI.199973     1    JLS -1.3643968

one_method_one_repeat_PGS_to_FDR_at_2 <- function(raw_PRS_onesimulation, FDR_threshold = 2){
  normalize_PRS <- function(x){
    x <- (x - mean(x))/ (sd(x))
    return(x)
  }
  
  raw_PRS_onesimulation[, PHENO:= factor(PHENO)]
  raw_PRS_onesimulation[, is.case := as.numeric((PHENO == 2))]
  raw_PRS_onesimulation[, pgs_score := normalize_PRS(pgs_score)]    
  
  FDR_by_threshold <- data.table()
  
  for(threshold in seq(-4, 4, length = 50)){
    raw_PRS_onesimulation[, predicted_case := as.numeric(pgs_score >= threshold)]
    temp_table <- raw_PRS_onesimulation[, .(true_positive = sum(predicted_case & is.case),
                                          false_positive = sum(predicted_case & !is.case),
                                          true_negative = sum(!predicted_case & !is.case),
                                          false_negative = sum(!predicted_case & is.case),
                                          threshold = threshold)]
    FDR_by_threshold <- rbind(FDR_by_threshold, temp_table)
  }
  
  FDR_by_threshold[, FDR := false_positive/(true_positive + false_positive)]
  FDR_by_threshold[, FPR:= false_positive/(false_positive + true_negative)]
  FDR_by_threshold[, FNR:= false_negative/(true_positive + false_negative)]
  FDR_by_threshold[, method:= raw_PRS_onesimulation$method[1]]
      
  FDR_go_to_table <- FDR_by_threshold[abs(threshold - FDR_threshold) < 0.05, ]
  return(FDR_go_to_table)
}


look_at_the_FDR_column <- one_method_one_repeat_PGS_to_FDR_at_2(raw_PRS_onesimulation, FDR_threshold = 2)

#THEN YOU CAN LOOP OVER ALL THE REPEATS FOR ONE METHOD AND STACK THOSE look_at_the_FDR_column TOGETHER
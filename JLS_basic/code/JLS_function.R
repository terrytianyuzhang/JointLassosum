
JLS_train <- function(JLS_population_weight_one,
                      JLS_l1_penalty,
                      JLS_shrinkage,
                      large_population_GWAS_file,
                      small_population_GWAS_file,
                      large_population_reference_prefix,
                      small_population_reference_prefix,
                      JLS_result_prefix,
                      large_population_type,
                      small_population_type,
                      chromosome = 1:22,
                      memory_limit = 2*10e9,
                      num_parallel_chr = 4){
  
  ####LOAD THE LD BLOCK BOUNDARY INFORMATION
  map_population_type_to_Berisa_label <- function(population_type){
    if(population_type %in% c("CEU", "EUR")){
      Berisa_label <- "EUR"
    }else if(population_type %in% c("YRI", "AFR")){
      Berisa_label <- "AFR"
    }else if(population_type %in% c("ASN")){
      Berisa_label <- "ASN"
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
  print(paste0('-----start model fitting for weight parameter = ', JLS_population_weight_one,'-----'))
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
                                mem.limit = memory_limit,
                                mc.cores = num_parallel_chr,
                                mc.preschedule = F)
  
  JLS_result_one_weight <- merge.mylassosum(JLS_result_by_chr)
  print(paste0('-----model fitting for weight parameter = ', JLS_population_weight_one,', finished-----'))
  
  ###OUTPUT LOCATION
  JLS_result_one_weight_file <- paste0(JLS_result_prefix, 
                                       sprintf("%.2f",JLS_population_weight_one), '.Rdata')
  save(JLS_result_one_weight,
       file=JLS_result_one_weight_file)
  
  JLS_regression_coefficient <- as.data.frame(JLS_result_one_weight$beta)
  names(JLS_regression_coefficient) <- paste0('l1_para_is_', JLS_l1_penalty)
  JLS_regression_coefficient$CHR <- pheno_gene_correlation$CHR
  JLS_regression_coefficient$ID <- pheno_gene_correlation$ID
  
  
  JLS_regression_coefficient_file <- paste0(JLS_result_prefix, 
                                            sprintf("%.2f",JLS_population_weight_one), '_coefficient.txt')
  fwrite(JLS_regression_coefficient, file = JLS_regression_coefficient_file)
}

split_chromosome_population_n_send_to_worker <- function(chromosome,
                                                         JLS_population_weight_one,
                                                         JLS_l1_penalty,
                                                         JLS_shrinkage,
                                                         large_population_type,
                                                         small_population_type,
                                                         LD_block_boundary,
                                                         reference_file,
                                                         pheno_gene_correlation,
                                                         mem.limit){
  
  pheno_gene_correlation_one_chr_large <- pheno_gene_correlation[pheno_gene_correlation$CHR == chromosome, large_population_type]
  pheno_gene_correlation_one_chr_small <- pheno_gene_correlation[pheno_gene_correlation$CHR == chromosome, small_population_type]
  
  rm(pheno_gene_correlation) #save memory
  
  LD_block_boundary_one_chr_large <- LD_block_boundary[[large_population_type]][LD_block_boundary[[large_population_type]]$chr == paste0("chr",chromosome),]
  LD_block_boundary_one_chr_large <- as.data.frame(LD_block_boundary_one_chr_large)
  LD_block_boundary_one_chr_small <- LD_block_boundary[[small_population_type]][LD_block_boundary[[small_population_type]]$chr == paste0("chr",chromosome),]
  LD_block_boundary_one_chr_small <- as.data.frame(LD_block_boundary_one_chr_small)
  
  JLS_result_one_chr <- mylassosum(cor1 = pheno_gene_correlation_one_chr_large, 
                                   cor2 = pheno_gene_correlation_one_chr_small,
                                   gamma = JLS_population_weight_one, 
                                   lambda = JLS_l1_penalty, 
                                   shrink = JLS_shrinkage,
                                   fileName1 = reference_file[[chromosome]][[large_population_type]], 
                                   fileName2 = reference_file[[chromosome]][[small_population_type]],
                                   LDblocks1=LD_block_boundary_one_chr_large, 
                                   LDblocks2=LD_block_boundary_one_chr_small,
                                   chunk = TRUE, 
                                   mem.limit = mem.limit,
                                   trace=2)
  # 
  # re=mylassosum(cor1=COR[COR$CHR == chr,gwasANC[1]], 
  #               cor2=COR[COR$CHR == chr,gwasANC[2]],
  #               fileName1 = referenceFiles[[chr]][[gwasANC[1]]], 
  #               fileName2 = referenceFiles[[chr]][[gwasANC[2]]],  
  #               gamma = gamma, lambda = lambda, shrink=shrink,
  #               chunk=TRUE, mem.limit=mem.limit,
  #               trace=2,
  #               LDblocks1=as.data.frame(LDblocks[[gwasANC[1]]][LDblocks[[gwasANC[1]]]$chr == paste0("chr",chr),]), 
  #               LDblocks2=as.data.frame(LDblocks[[gwasANC[2]]][LDblocks[[gwasANC[2]]]$chr == paste0("chr",chr),]))
  return(JLS_result_one_chr)
}

mylassosum <- function(cor1, cor2,  fileName1, fileName2,
                       gamma=0.5,
                       lambda=exp(seq(log(0.001), log(0.05), length.out=5)), 
                       shrink=0.9, 
                       #ridge=F, thr= 10^(-4), 
                       chunk=T,
                       init=NULL, trace=0, maxiter=10000,
                       keep1=NULL, keep2=NULL, remove=NULL, extract1=NULL, extract2=NULL, exclude=NULL, distance=NULL,
                       chr=NULL, mem.limit=4*10^9,
                       LDblocks1=NULL, LDblocks2=NULL) {
  
  cor1 <- as.vector(cor1)
  cor2 <- as.vector(cor2)
  stopifnot(!any(is.na(cor1)))
  stopifnot(!any(is.na(cor2)))
  stopifnot(length(cor1) == length(cor2))
  
  print(paste0("distance is NULL:", is.null(distance)))
  parsed1<-parseselect(fileName1, extract=extract1, exclude=exclude, keep=keep1, remove=remove, distance=distance, chr=chr) #extract1 could be a logical vector? 
  parsed2<-parseselect(fileName2, extract=extract2, exclude=exclude, keep=keep2, remove=remove, distance=distance, chr=chr)
  
  ### Split/group into chunks###
  if(chunk==T){
    required.memory=as.numeric(parsed1$p) * as.numeric(parsed1$n) * 8 + as.numeric(parsed2$p) * as.numeric(parsed2$n) *8
    n.chunks<-ceiling(required.memory/mem.limit)
    print(paste0("Number of chunks: ",n.chunks))
    
    chunk.len1=ceiling(parsed1$p/n.chunks)
    chunk.len2=ceiling(parsed2$p/n.chunks)
    stopifnot(chunk.len1 == chunk.len2)
    chunk.len=chunk.len1
    print(paste0("chunk length: ",chunk.len))
    
    if(n.chunks > 1) {
      if(is.null(parsed1$extract)){
        extract1=lapply(1:n.chunks, function(i) logical.vector(min(((i-1)*chunk.len + 1), parsed1$P) : min((i*chunk.len), parsed1$P), parsed1$P))
      } else{
        extract1=lapply(1:n.chunks, function(i) logical.vector(parsed1$extract & (cumsum(parsed1$extract)<=i*chunk.len) & (cumsum(parsed1$extract) > (i-1)*chunk.len),parsed1$P))
      }
      if(is.null(parsed2$extract)){
        extract2=lapply(1:n.chunks, function(i) logical.vector(min(((i-1)*chunk.len + 1), parsed2$P) : min((i*chunk.len), parsed2$P), parsed2$P))
      } else{
        extract2=lapply(1:n.chunks, function(i) logical.vector(parsed2$extract & (cumsum(parsed2$extract)<=i*chunk.len) & (cumsum(parsed2$extract) > (i-1)*chunk.len),parsed2$P))
      }
      
      results.list <- lapply(1:n.chunks, function(i) {
        mylassosum(cor1=cor1, cor2=cor2, fileName1=fileName1, fileName2=fileName2, 
                   gamma=gamma, lambda=lambda, shrink=shrink, 
                   keep1=parsed1$keep, keep2=parsed2$keep, extract1=extract1[[i]],extract2=extract2[[i]], distance = NULL,
                   trace=trace, maxiter=maxiter, chunk=FALSE,
                   LDblocks1=LDblocks1, LDblocks2=LDblocks2,  
                   mem.limit=mem.limit)
      })
      print(paste0("length(results.list): ", length(results.list)) )
      print(paste0("class: ", sapply(results.list, "class") ))
      
      re.lassosum=merge.mylassosum(results.list)
      #      return(do.call("merge.mylassosum", results.list))
      return(re.lassosum)
    }
  }
  
  if(!is.null(parsed1$extract) & length(parsed1$extract) == length(cor1)) cor1=cor1[parsed1$extract]
  if(!is.null(parsed2$extract) & length(parsed2$extract) == length(cor2)) cor2=cor2[parsed2$extract]
  ### Split by LD region ###
  reindex <- function(x) {
    as.integer(x)-min(as.integer(x)) + 1
  }
  
  bimfile1=read.table2(paste0(fileName1,".bim"))
  bimfile2=read.table2(paste0(fileName2,".bim"))
  
  bimfile1$V1 = paste0("chr",bimfile1$V1)
  bimfile2$V1 = paste0("chr",bimfile2$V1)
  
  #  if(is.null(extract)){
  #    INDEX1 = factor(bimfile1$V1, levels=unique(bimfile1$V1))
  #    INDEX2 = factor(bimfile2$V1, levels=unique(bimfile2$V1))
  #  } else {
  #    INDEX1 = factor(bimfile1$V1[parsed1$extract], levels=unique(bimfile1$V1[parsed1$extract]))
  #    INDEX2 = factor(bimfile2$V1[parsed2$extract], levels=unique(bimfile2$V1[parsed2$extract]))
  #  }
  
  
  if(is.null(extract1)){
    blocks1 <- splitgenome(CHR = bimfile1$V1, 
                           POS = bimfile1$V4,
                           ref.CHR = LDblocks1[,1], 
                           ref.breaks = LDblocks1[,3])
  } else {
    blocks1 <- splitgenome(CHR = bimfile1$V1[parsed1$extract], 
                           POS = bimfile1$V4[parsed1$extract],
                           ref.CHR = LDblocks1[,1], 
                           ref.breaks = LDblocks1[,3])
  }
  
  if(is.null(extract2)){
    blocks2 <- splitgenome(CHR = bimfile1$V1, 
                           POS = bimfile1$V4,
                           ref.CHR = LDblocks2[,1], 
                           ref.breaks = LDblocks2[,3])
  } else {
    blocks2 <- splitgenome(CHR = bimfile1$V1[parsed2$extract], 
                           POS = bimfile1$V4[parsed2$extract],
                           ref.CHR = LDblocks2[,1], 
                           ref.breaks = LDblocks2[,3])
  }
  
  #blocks1<-as.numeric(unlist(tapply(blocks1, INDEX1, reindex)))
  #blocks1<-as.integer(blocks1)-min(as.integer(blocks1)) + 1
  blocks1<-as.numeric(factor(blocks1, levels=unique(blocks1)))
  blocks2<-as.numeric(factor(blocks2, levels=unique(blocks2)))
  #blocks2<-as.integer(blocks2)-min(as.integer(blocks2)) + 1
  #blocks2<-as.numeric(unlist(tapply(blocks2, INDEX2, reindex)))
  
  el <- myelnet(lambda, shrink, gamma, fileName1, fileName2, cor1, cor2, 
                #thr=thr,
                parsed1, parsed2,
                trace=trace, maxiter=maxiter, 
                blocks1=blocks1, blocks2=blocks2)
  
  #--------------------------------------------
  #This section is not used yet
  #if(ridge) {
  #  if(is.null(blocks)) blocks <- rep(1, p)
  #  reps <- 1:max(blocks)
  #  svd <- lapply(1:reps, function(i) svd(X[,blocks==i], nu=0))
  #  invert <- lapply(1:reps, function(i) 1/(svd[[i]]$d^2 + shrink))
  #  VG <- lapply(1:reps, function(i) svd[[i]]$v %*% Diagonal(x=invert[[i]] - 1/shrink))
  #  VTr <- lapply(1:reps, function(i) t(svd[[i]]$v) %*% cor[blocks == i])
  #  VGVTr <- lapply(1:reps, function(i) VG[[i]] %*% VTr[[i]])
  #  Ridge <- lapply(1:reps, function(i) as.vector(VGVTr[[i]] + 1/shrink * cor[blocks==i]))
  #  Ridge <- unlist(Ridge)
  #} else Ridge <- NULL
  #----------------------------------------------
  
  nparams <- colSums(el$beta != 0) 
  
  toreturn <- list(lambda=lambda, 
                   beta=el$beta,
                   conv=el$conv,
                   pred1=el$pred1,
                   pred2=el$pred2,
                   loss=el$loss,
                   trainerror1 = el$trainerror1,
                   trainerror2 = el$trainerror2,
                   fbeta=el$fbeta,
                   sd1=el$sd1,
                   sd2=el$sd2,
                   shrink=shrink,
                   nparams=nparams, 
                   gamma=gamma
                   #,ridge=Ridge
  )
  
  class(toreturn) <- "lassosum"
  
  
  return(toreturn)
  
  #' @return A list with the following
  #' \item{lambda}{same as the lambda input}
  #' \item{beta}{A matrix of estimated coefficients}
  #' \item{conv}{A vector of convergence indicators. 1 means converged. 0 not converged.}
  #' \item{pred}{\eqn{=(1-s)X\beta}}
  #' \item{loss}{\eqn{=(1-s)\beta'X'X\beta/n - 2\beta'r}}
  #' \item{trainerror1}{training error for population 1}
  #' \item{trainerror2}{training error for population 2}
  #' \item{fbeta}{\eqn{=\beta'R\beta - 2\beta'r + 2\lambda||\beta||_1}}
  #' \item{sd}{The standard deviation of the reference panel SNPs}
  #' \item{shrink}{same as input}
  #' \item{nparams}{Number of non-zero coefficients}
  #' \item{gamma}{A parameter balancing two populations}
  #' \item{ridge}{ridge regression estimates}
  
}

#' @title Elastic net using summary statistics
#' @description Coordinate descent algorithm to solve: 
#' 0.5 x'X'Xx - x'b + lambda1 ||x||_1 + 0.5 lambda2 ||x||_2^2
#' Function to get elastic net solutions given X, a reference panel, and
#' b, regression coefficients
#' @keywords internal
myelnet <- function(lambda1, lambda2=0.0, gamma, fileName1, fileName2, b1, b2, 
                    #thr=1e-4, #lambda2 is s
                    parsed1, parsed2,
                    trace=0, maxiter=10000, 
                    blocks1=NULL, blocks2=NULL,
                    x=NULL) {
  #stopifnot(length(b1) == ncol(X1))
  #stopifnot(length(b2) == ncol(X2))
  #stopifnot(length(b1) == length(b2))
  #diag1 <- colSums(X1^2)
  #diag2 <- colSums(X2^2)
  
  N1 <- parsed1$N
  N2 <- parsed2$N
  stopifnot(parsed1$p==parsed2$p)
  P <- parsed1$P
  
  if(length(lambda2) > 1) {
    nlambda2 <- length(lambda2)
    for(i in 1:nlambda2) {
      result <- myelnet(lambda1, lambda2[i], gamma, fileName1, fileName2, b1, b2, 
                        #thr,
                        parsed1, parsed2,
                        trace, maxiter, x)
      result <- list(fit=result, lambda2=lambda2[i])
      if(i == 1) Result <- rep(result, nlambda2) else
        Result[i] <- result
      
    }
    return(Result)
  }
  
  #order <- order(lambda1, decreasing = T)
  #lambda1a <- lambda1[order]
  #conv <- lambda1a * NA
  len <- length(b1)
  #beta <- matrix(NA, len, length(lambda1))
  #pred1 <- matrix(NA, N1, length(lambda1))
  #pred2 <- matrix(NA, N2, length(lambda1))
  #loss <- rep(NA, length(lambda1))
  #fbeta <- loss
  
  if(is.null(x)) x <- b1 * 0.0 else { #initialize beta with 0s
    stopifnot(length(x) == len)
    x <- x + 0.0 # Making sure R creates a copy...
  }
  
  if(is.null(blocks1)) {
    Blocks1 <- list(startvec=0, endvec=len - 1)
    blocks1 <- rep(0, len)
  } else {
    Blocks1 <- parseblocks(blocks1)
    blocks1 <- blocks1 - 1 # for in C, pos starts from 0
    stopifnot(max(Blocks1$endvec)==len - 1) # leave it hear, the functions should be the same
  }
  
  if(is.null(blocks2)) {
    Blocks2 <- list(startvec=0, endvec=len - 1)
    blocks2 <- rep(0, len)
  } else {
    Blocks2 <- parseblocks(blocks2)
    blocks2 <- blocks2 - 1 # for in C, pos starts from 0
    stopifnot(max(Blocks2$endvec)==len - 1) # leave it hear, the functions should be the same
  }
  
  if(is.null(parsed1$extract)) {
    extract1 <- list(integer(0), integer(0))
    #col_skip_pos = integer(0); col_skip = integer(0)
  } else {
    #print(parsed1$extract)
    extract1 <- selectregion(!parsed1$extract)
    extract1[[1]] <- extract1[[1]] - 1
  }
  
  if(is.null(parsed1$keep)) {
    keepbytes1 <- integer(0)
    keepoffset1 <- integer(0)
  } else {
    pos <- which(parsed1$keep) - 1
    keepbytes1 <- floor(pos/4)
    keepoffset1 <- pos %% 4 * 2
  }
  
  if(is.null(parsed2$extract)) {
    extract2 <- list(integer(0), integer(0))
    #col_skip_pos = integer(0); col_skip = integer(0)
  } else {
    #print(parsed2$extract)
    extract2 <- selectregion(!parsed2$extract)
    extract2[[1]] <- extract2[[1]] - 1
  }
  
  if(is.null(parsed2$keep)) {
    keepbytes2 <- integer(0)
    keepoffset2 <- integer(0)
  } else {
    pos <- which(parsed2$keep) - 1
    keepbytes2 <- floor(pos/4)
    keepoffset2 <- pos %% 4 * 2
  }
  #X1 <- as.matrix(X1)
  #X2 <- as.matrix(X2)
  #yhat1 <- as.vector(X1 %*% x) #initialize PRS pred for population 1
  #yhat2 <- as.vector(X2 %*% x) #initialize PRS pred for population 2
  
  
  #  for(i in 1:length(lambda1a)) {
  #    if(trace > 0) cat("lambda1: ", lambda1a[i], "\n")
  #conv[i] <- myrepelnet(lambda1a[i], lambda2, gamma, diag1, diag2, X1, X2, b1, b2, 
  #thr,
  #                      x,yhat1, yhat2, trace-1,maxiter,
  #                      Blocks1$startvec, Blocks1$endvec, Blocks2$startvec, Blocks2$endvec, blocks1, blocks2)
  #if(conv[i] != 1) stop("Not converging...")
  
  ####This is a C worker function
  results<- myrunElnet(lambda=lambda1, shrink=lambda2, gamma=gamma, 
                       fileName1=fileName1, fileName2=fileName2,
                       r1=b1, r2=b2, 
                       N1=N1, N2=N2, P=P, 
                       col_skip_pos1 = extract1[[1]], col_skip1 = extract1[[2]], 
                       keepbytes1 = keepbytes1, keepoffset1 = keepoffset1, 
                       col_skip_pos2 = extract2[[1]], col_skip2 = extract2[[2]], 
                       keepbytes2 = keepbytes2, keepoffset2 = keepoffset2, 
                       #thr, 
                       x, trace, maxiter, 
                       startvec1 = Blocks1$startvec, endvec1 = Blocks1$endvec,
                       startvec2 = Blocks2$startvec, endvec2 = Blocks2$endvec,
                       blocks1, blocks2)
  
  #    beta[,i] <- x
  #    pred1[,i] <- yhat1
  #    pred2[,i] <- yhat2
  #    loss[i] <- gamma * (sum(yhat1^2) - 2* sum(b1 * x)) + (1-gamma) * (sum(yhat2^2) - 2* sum(b2 * x)) 
  #    fbeta[i] <- loss[i] + 2* sum(abs(x))*lambda1a[i] + sum(x^2)*lambda2
  #  }
  
  
  #conv[order] <- conv
  #beta[,order] <- beta
  #pred1[,order] <- pred1
  #pred2[,order] <- pred2
  #loss[order] <- loss
  #fbeta[order] <- fbeta
  
  return(list(lambda1=lambda1, lambda2=lambda2, gamma=gamma, beta=results$beta, conv=results$conv, 
              pred1=results$pred1, pred2=results$pred2, loss=results$loss, fbeta=results$fbeta,
              trainerror1 = results$trainerror1, trainerror2 = results$trainerror2, 
              sd1=results$sd1, sd2=results$sd2))
  
}

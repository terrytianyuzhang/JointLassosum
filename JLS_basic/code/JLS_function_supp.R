#' @title Merge lassosum results 
#' @description e.g. when calculated over different blocks/chromosomes
#' @method merge lassosum
#' @export
#' 
merge.mylassosum <- function(ll) {
  
  #  ll <- list(...)
  #  stopifnot(all(sapply(ll, "class") == "lassosum"))
  stopifnot(all(sapply(ll, function(x) all(x$lambda == ll[[1]]$lambda))))
  shrink <- sapply(ll, function(x) x$shrink)
  print('shrink is')
  print(shrink)
  print('shrink[1] is')
  print(shrink[1])
  print('all(shrink == shrink[1])')
  print(shrink == shrink[1])
  print(all(shrink == shrink[1]))
  stopifnot(all(shrink == shrink[1]))
  gamma <- sapply(ll, function(x) x$gamma)
  stopifnot(all(gamma == gamma[1]))
  
  Cumsum <- function(...) {
    mat <- do.call("rbind", list(...))
    if(ncol(mat) > 0) return(as.vector(colSums(mat))) else
      return(numeric(0))
  }
  results <- list()
  results$lambda <- ll[[1]]$lambda
  results$beta <- do.call("rbind", lapply(ll, function(x) x$beta))
  results$conv <- do.call("pmin", lapply(ll, function(x) x$conv))
  pred1 <- do.call("Cumsum", lapply(ll, function(x) as.vector(x$pred1)))
  pred2 <- do.call("Cumsum", lapply(ll, function(x) as.vector(x$pred2)))
  results$pred1 <- matrix(pred1, ncol=length(results$lambda), nrow=nrow(ll[[1]]$pred1))
  results$pred2 <- matrix(pred2, ncol=length(results$lambda), nrow=nrow(ll[[1]]$pred2))
  results$loss <- do.call("Cumsum", lapply(ll, function(x) x$loss))
  results$trainerror1 <- do.call("Cumsum", lapply(ll, function(x) matrix(x$trainerror1, nrow = 1)))
  results$trainerror2 <- do.call("Cumsum", lapply(ll, function(x) matrix(x$trainerror2, nrow = 1)))
  
  results$fbeta <- do.call("Cumsum", lapply(ll, function(x) x$fbeta))
  results$sd1 <- do.call("c", lapply(ll, function(x) x$sd1))
  results$sd2 <- do.call("c", lapply(ll, function(x) x$sd2))
  results$shrink <- ll[[1]]$shrink
  results$nparams <- do.call("Cumsum", lapply(ll, function(x) x$nparams))
  results$gamma <- ll[[1]]$gamma
  class(results) <- "lassosum"
  return(results)
}

#' @title Obtains the number of column (SNPs) in a PLINK bfile
#' 
#' @param bfile Plink file stem
#' @return an integer with the number of columns
#' #@keywords internal
#' @export
ncol.bfile <- function(bfile) {
  bimfile <- paste0(bfile, ".bim")
  if(!file.exists(bimfile)) 
    stop(paste0("Cannot find ", bimfile)) 
  
  return(countLines(bimfile))
  # if(Sys.info()["sysname"] == "Windows") {
  # 	wc.output <- shell(paste("wc -l", bimfile), intern=T)
  # } else {
  # 	wc.output <- system(paste("wc -l", bimfile), intern=T)
  # }
  # return(as.numeric(strsplit(wc.output, split = "\\s+")[[1]][1]))
}

#' @title Obtains the number of individuals in a PLINK bfile
#' 
#' @param bfile Plink file stem
#' @return an integer with the number of rows
#' #@keywords internal
#' @export
nrow.bfile <- function(bfile) {
  famfile <- paste0(bfile, ".fam")
  if(!file.exists(famfile)) 
    stop(paste0("Cannot find ", famfile)) 
  return(countLines(famfile))
  # if(Sys.info()["sysname"] == "Windows") {
  # 	wc.output <- shell(paste("wc -l", famfile), intern=T)
  # } else {
  # 	wc.output <- system(paste("wc -l", famfile), intern=T)
  # }
  # return(as.numeric(strsplit(wc.output, split = "\\s+")[[1]][1]))
}

parse.pheno.covar <- function(pheno, covar, parsed, trace=0) {
  #' @keywords internal
  fam <- parsed[['fam']]
  keep <- parsed$keep
  # keep <- NULL
  pheno.df <- NULL
  
  update.keep <- function(old, new) {
    if(all(new)) {
      return(old)
    } else {
      if(is.null(old)) return(new) else {
        if(is.null(new)) return(old) else 
          return(old & new)
      }
    }
  }
  #### pheno ####
  if(!is.null(pheno) && is.character(pheno) && length(pheno) == 1) {
    if(file.exists(pheno)) pheno <- read.table2(pheno, header=T) else 
      stop(paste("Cannot find", pheno))
  }
  if(is.data.frame(pheno)) {
    if(ncol(pheno) != 3) {
      stop(paste("A pheno data.frame must have 3 columns exactly",
                 "with the first 2 with headers 'FID' and 'IID'"))
    }
    colnames <- colnames(pheno) 
    if(!all(colnames[1:2] == c("FID", "IID"))) {
      stop(paste("The first two columns of the pheno", 
                 "data.frame must have headers 'FID' and 'IID'"))
    }
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    pheno.df <- pheno
    colnames(pheno.df)[3] <- "pheno"
    rownames(pheno) <- paste(pheno$FID, pheno$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(pheno))
    Pheno <- as.data.frame(pheno)[,3]
    names(Pheno) <- rownames(pheno)
  } else {
    if(!is.null(pheno)) {
      stopifnot(length(pheno) == parsed$n)
    } else {
      fam <- read.table2(parsed$famfile)
      if(is.null(parsed$keep)) pheno <- fam$V6 else 
        pheno <- fam$V6[parsed$keep]
    }
  }
  
  #### covar ####
  user.covar <- FALSE
  if(!is.null(covar) && is.character(covar) && length(covar) == 1) {
    if(file.exists(covar)) covar <- read.table2(covar, header=T) else 
      stop(paste("Cannot find", covar))
  }
  if(is.data.frame(covar) & all(colnames(covar)[1:2] == c("FID", "IID"))) {
    user.covar <- TRUE
    covar <- as.data.frame(covar)
    colnames <- colnames(covar) 
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    rownames(covar) <- paste(covar$FID, covar$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(covar))
    Covar <- covar[,-(1:2), drop=FALSE]
  } else {
    if(!is.null(covar)) {
      if(is.vector(covar)) covar <- matrix(covar, ncol=1)
      if(is.matrix(covar)) covar <- as.data.frame(covar)
      Covar <- covar
    } 
  }
  
  #### updates ####
  parsed$keep <- update.keep(parsed$keep, keep)
  if(!is.null(parsed$keep)) parsed$n <- sum(parsed$keep)
  if(is.data.frame(pheno)) {
    if(!is.null(parsed$keep)) {
      names <- rownames(fam)[parsed$keep]
    } else {
      names <- rownames(fam)
    }
    pheno <- Pheno[names] 
    if(trace) {
      message(length(pheno), " out of ", length(Pheno), " samples kept in pheno.")
      # message(paste("Note that the order of best.pgs is the order given in the .fam file", 
      #               " rather than the pheno data.frame. Use v$best.pgs[v$order] to get", 
      #               " the pgs in the order of the phenotype."))
    }
    Order <- 1:length(pheno)
    names(Order) <- names
    pheno.df$order <- Order[names(Pheno)]
  } 
  
  if(user.covar) {
    if(!is.null(parsed$keep)) covar <- Covar[rownames(fam)[parsed$keep],,drop=F] else 
      covar <- Covar[rownames(fam),,drop=F]
    if(trace) message(nrow(covar), " out of ", nrow(Covar), " samples kept in covar.")
  } 
  
  if(length(pheno) == 0) {
    stop("No phenotype left. Perhaps the FID/IID do not match?")
  } else if(length(pheno) != parsed$n) {
    stop("The length of pheno does not match the number of samples.")
  }
  if(!is.null(covar) && nrow(covar) != parsed$n) {
    stop(paste("The dimension of the covar matrix does not match the number of samples used.", 
               "If your covariate is a data.frame with FID and IID, make sure they have headers."))
  }
  # if(sd(pheno, na.rm = TRUE) == 0) stop("There's no variation in phenotype")
  parsed$fam <- fam
  
  return(list(pheno=pheno, covar=covar, parsed=parsed, table=pheno.df))
  
}

parseblocks <- function(vec) {
  
  #' @keywords internal
  if(is.factor(vec)) vec <- as.integer(vec)
  vec <- as.integer(factor(vec, levels=unique(vec)))
  blocks <- unique(vec)
  stopifnot(blocks == sort(blocks))
  stopifnot(min(blocks) == 1)
  stopifnot(max(blocks) == length(blocks))
  rle <- rle(vec)
  endvec <- cumsum(rle$lengths)
  startvec <- c(0, endvec[-length(endvec)])
  endvec <- endvec - 1 
  return(list("startvec"=startvec, "endvec"=endvec))
  
}

#' @title Parse the keep/remove/extract/exclude/chr options
#' @details \code{keep} and \code{remove} can take one of three
#' formats:
#' 1. A logical vector indicating which indivduals to keep/remove,
#'
#' 2. A \code{data.frame} with two columns giving the FID and IID of the indivdiuals to keep/remove (matching those in the .fam file), or
#' 3. A character scalar giving the text file with the FID/IID.
#'
#' Note that these files should not contain headers.
#'
#' \code{extract} and \code{exclude} can be of format type __(1)__ or __(3)__ describe above, or a character vector of SNP ids (matching those in the .bim file).

#' @param bfile plink file stem
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param export Logical. Whether to include the \code{bim} and/or \code{fam} data in the returned object.
#' @param order.important Logical. Whether to strictly enforce that the order of \code{extract} and/or \code{keep} must match with \code{bim} or \code{fam}, respectively.
#'
#' @md
#' @export
parseselect <- function(bfile, extract=NULL, exclude=NULL,
                        keep=NULL, remove=NULL, chr=NULL,
                        export=FALSE, distance=NULL, order.important=FALSE) {
  
  #### Introduced for parsing multiple bfiles ####
  if(length(bfile) > 1) {
    if(!is.null(chr) || !is.null(extract) || !is.null(exclude)) {
      stop("bfile cannot be a vector if we are parsing extract/exclude/chr. yet.")
    }
    bfile <- bfile[1]
  }
  
  #### Checks ####
  stopifnot(is.character(bfile) && length(bfile) == 1)
  bedfile <- paste0(bfile, ".bed")
  bimfile <- paste0(bfile, ".bim")
  famfile <- paste0(bfile, ".fam")
  stopifnot(file.exists(bedfile))
  stopifnot(file.exists(bimfile))
  stopifnot(file.exists(famfile))
  
  if(grepl("^~", bfile)) {
    stop("Don't use '~' as a shortcut for the home directory.")
  }
  
  p <- P <- ncol.bfile(bfile)
  n <- N <- nrow.bfile(bfile)
  bim <- NULL
  fam <- NULL
  
  #### extract ####
  if(!is.null(extract)) {
    if(is.logical(extract)) {
      stopifnot(length(extract) == P)
    } else {
      if(is.null(attr(extract, "not.a.file", exact=T))) {
        if(is.character(extract) && length(extract) == 1 && file.exists(extract)) {
          ### I'm interpreting this as a filename
          SNPs <- read.table2(extract)
          stopifnot(ncol(SNPs)==1)
          extract <- SNPs[[1]]
        }
      } else {
        attr(extract, "not.a.file") <- NULL
      }
      if(is.vector(extract)) {
        if(is.null(distance)){
          Extract <- as.character(extract)
          bim <- read.table2(bimfile)
          extract <- bim$V2 %in% Extract
          if(order.important) {
            if(!all(bim$V2[extract] == Extract)) {
              stop("Order of extract SNPs does not match that in .bim file.")
            }
          }
        } else if(distance>0){
          Extract <- as.character(extract)
          bim <- read.table2(bimfile)
          bim.extract <- bim[bim$V2 %in% Extract,]
          bim.extract$start = bim.extract$V4 - distance
          bim.extract$start = ifelse(bim.extract$start<0, 0, bim.extract$start)
          bim.extract$end = bim.extract$V4 + distance
          extract<-rep(FALSE,nrow(bim))
          for(l in 1:nrow(bim.extract)){
            extract[bim$V4 > bim.extract$start[l] & bim$V4 < bim.extract$end[l] & bim.extract$V1[l]==bim$V1]<-TRUE
          }
        } else if(distance==0){
          Extract <- as.character(extract)
          bim <- read.table2(bimfile)
          extract <- bim$V2 %in% Extract
          if(order.important) {
            if(!all(bim$V2[extract] == Extract)) {
              stop("Order of extract SNPs does not match that in .bim file.")
            }
          }
        } else{
          stop("distance should be a non negative number ")
        }
      } else {
        stop("I don't know what to do with this type of input for extract")
      }
    }
    
    p <- sum(extract)
  }
  
  #### exclude ####
  if(!is.null(exclude)) {
    if(is.logical(exclude)) {
      stopifnot(length(exclude) == P)
    } else {
      if(is.null(attr(exclude, "not.a.file", exact=T))) {
        if(is.character(exclude) && length(exclude) == 1 && file.exists(exclude)) {
          ### I'm interpreting this as a filename
          SNPs <- read.table2(exclude)
          stopifnot(ncol(SNPs)==1)
          exclude <- SNPs[[1]]
        }
      } else {
        attr(exclude, "not.a.file") <- NULL
      }
      if(is.vector(exclude)) {
        exclude <- as.character(exclude)
        if(is.null(bim)) bim <- read.table2(bimfile)
        exclude <- bim$V2 %in% exclude
      } else {
        stop("I don't know what to do with this type of input for exclude")
      }
      
    }
    
    if(is.null(extract)) extract <- !exclude else
      extract <- extract & !exclude
    
    p <- sum(extract)
  }
  
  #### chr ####
  if(!is.null(chr)) {
    
    stopifnot(is.vector(chr))
    chr <- as.character(chr)
    
    if(is.null(bim)) bim <- read.table2(bimfile)
    bimchr <- bim$V1
    bimchr[bimchr==""]
    extract.chr <- bim$V1 %in% chr
    
    if(is.null(extract)) extract <- extract.chr else
      extract <- extract & extract.chr
    
    p <- sum(extract)
    
  }
  
  #### keep ####
  if(!is.null(keep)) {
    if(is.logical(keep)) {
      stopifnot(length(keep) == N)
    } else {
      if(is.character(keep) && length(keep) == 1 && file.exists(keep)) {
        ### I'm interpreting this as a filename
        keep <- read.table2(keep)
      }
      if(is.vector(keep)) {
        #keep <- as.data.frame(keep)
        #stopifnot(ncol(keep)==2)
        fam <- read.table2(famfile)
        famID <- fam[,1]
        #keepID <- keep
        keep <- famID %in% keep
        if(order.important) {
          if(!all(famID[keep] == keepID)) {
            stop("Order of keep doesn't match that in .fam file")
          }
        }
      } else {
        stop("I don't know what to do with this type of input for keep")
      }
      
    }
    n <- sum(keep)
  }
  
  #### remove ####
  if(!is.null(remove)) {
    if(is.logical(remove)) {
      stopifnot(length(remove) == N)
    } else {
      if(is.character(remove) && length(remove) == 1 && file.exists(remove)) {
        ### I'm interpreting this as a filename
        remove <- read.table2(remove)
      }
      if(inherits(remove, "data.frame")) {
        remove <- as.data.frame(remove)
        stopifnot(ncol(remove)==2)
        if(is.null(fam)) fam <- read.table2(famfile)
        famID <- paste(fam[,1], fam[,2], sep=".")
        removeID <- paste(remove[,1], remove[,2], sep=".")
        remove <- famID %in% removeID
      } else {
        stop("I don't know what to do with this type of input for remove")
      }
    }
    
    if(is.null(keep)) keep <- !remove else
      keep <- keep & !remove
    
    n <- sum(keep)
  }
  
  if(n==0) stop("No individuals left after keep/remove! Make sure the FID/IID are correct.")
  if(p==0) stop("No SNPs left after extract/exclude/chr! Make sure the SNP ids are correct.")
  
  if(!export) {
    return(list(keep=keep, extract=extract,
                N=N, P=P, n=n, p=p, bfile=bfile,
                bimfile=bimfile, famfile=famfile,
                bim=NULL, fam=NULL))
  } else {
    return(list(keep=keep, extract=extract,
                N=N, P=P, n=n, p=p, bfile=bfile,
                bimfile=bimfile, famfile=famfile,
                bim=bim, fam=fam))
  }
  #' @return A list containing:
  #' \item{keep}{Either NULL or a logical vector of which individuals to keep}
  #' \item{extract}{Either NULL or a logical vector of which SNPs to extract}
  #' \item{N}{Number of rows in the PLINK bfile}
  #' \item{P}{Number of columns in the PLINK bfile}
  #' \item{n}{Number of rows in the PLINK bfile after keep}
  #' \item{p}{Number of columns in the PLINK bfile after extract}
  #' \item{bfile}{File path to bfile stub}
  #' \item{bimfile}{File path to bim file}
  #' \item{famfile}{File path to fam file}
  #' \item{bim}{Either \code{NULL} or a data frame of bim data}
  #' \item{fam}{Either \code{NULL} or a data frame of fam data}
  
}

#' @title Function to read a text file 
#' @keywords internal
read.table2 <- function(file, header=F, data.table=F, check.names=TRUE, ...) {
  return(data.table::fread(file, header=header, data.table=data.table, 
                           check.names=check.names, ...))
}

#' @title Internal function to parse extract
#'
#' @param keep A boolean vector of which position to keep
#' @keywords internal
#' 
selectregion <- function(keep) {
  
  pos.keep <- which(keep)
  end <- pos.keep[-1]
  start <- pos.keep[-length(pos.keep)]
  diff.keep <- end - start
  skip <- diff.keep > 1
  starts <- c(start[1], end[skip])
  rle <- rle(keep)
  # print(rle)
  lengths <- rle$lengths[rle$values]
  return(list(starts, lengths))
  
}

splitgenome <- function(CHR, POS, ref.CHR, ref.breaks, details=T, right=TRUE) {
  #' Function to split a set of SNPs by their position using a reference
  #' Break points should be given for the reference by chromosome
  #' @keywords internal
  
  CHR <- as.character(CHR)
  ref.CHR <- as.character(ref.CHR)
  POS <- as.integer(POS)
  ref.breaks <- as.integer(ref.breaks)
  
  stopifnot(all(!is.na(POS)) && all(!is.na(ref.breaks)))
  stopifnot(all(POS >= 1) && all(ref.breaks >= 1))
  stopifnot(length(CHR) == length(POS))
  stopifnot(length(ref.CHR) == length(ref.breaks))
  
  
  chr <- (unique(CHR))
  chr.ref <- (unique(ref.CHR))
  included.chr <- chr %in% chr.ref
  # if(all(!included.chr)) stop("Cannot match the chromosomes. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
  if(!all(included.chr)) stop("Some chromosomes were not defined in the reference. Make sure the notations are the same. e.g. 'chr1' vs 1 or chrX vs chr23.")
  
  levels <- character(0)
  results <- character(length(POS))
  Details <- data.frame()
  for(C in chr.ref) {
    breaks <- sort(unique(ref.breaks[ref.CHR == C]))
    if(breaks[1] > 1) breaks <- c(1, breaks)
    if(breaks[length(breaks)] < Inf) breaks <- c(breaks, Inf)
    cut <- cut(POS[CHR == C],include.lowest = T,breaks = breaks, right=right)
    levels <- c(levels, paste0(C, "_", levels(cut)))
    cut <- paste0(C, "_", as.character(cut))
    results[CHR == C] <- cut
    if(details) {
      df <- data.frame(chr=C, start=breaks[-length(breaks)], end=breaks[-1])
      Details <- rbind(Details, df)
    }
  }
  results <- factor(results, levels = levels)
  
  if(details) {
    Details$counts <- as.integer(table(results))
    attr(results, "details") <- Details
  }
  
  return(results)
  
}

### Validation function ###
validate<-function(fit, fileName, extract=NULL, keep=NULL, distance=NULL, chunk=T, mem.limit=4*10^9,...){
  
  parsed<-parseselect(fileName, extract=extract, keep=keep, distance=distance)
  ### Split/group into chunks###
  if(chunk==T){
    required.memory=as.numeric(parsed$p) * as.numeric(parsed$n) * 8 
    n.chunks<-ceiling(required.memory/mem.limit)
    print(paste0("Number of chunks: ",n.chunks))
    
    chunk.len=ceiling(parsed$p/n.chunks)
    print(paste0("chunk length: ",chunk.len))
    
    if(n.chunks > 1) {
      if(is.null(parsed$extract)){
        extract=lapply(1:n.chunks, function(i) logical.vector(min(((i-1)*chunk.len + 1), parsed$P) : min((i*chunk.len), parsed$P), parsed$P))
        fit.list<-vector("list",length=n.chunks)
        for(i in 1:n.chunks) {
          fit.list[[i]] = fit
          fit.list[[i]]$beta=fit.list[[i]]$beta[extract[[i]],,drop=FALSE]
        }
      } else{
        extract=lapply(1:n.chunks, function(i) logical.vector(which((cumsum(parsed$extract)<=i*chunk.len) & (cumsum(parsed$extract) > (i-1)*chunk.len)),parsed$P))
        fit.list<-vector("list",length=n.chunks)
        for(i in 1:n.chunks) {
          fit.list[[i]] = fit
          fit.list[[i]]$beta = fit.list[[i]]$beta[which(parsed$extract) %in% which(extract[[i]]),,drop=FALSE]
        }
      }
      
      results.list <- lapply(1:n.chunks, function(i) {
        print(i)
        validate(fit=fit.list[[i]], fileName = fileName, extract=extract[[i]], chunk = FALSE, distance=NULL)
      })
      
      mysum <- function(...) {
        mat <- do.call("rbind", list(...))
        if(ncol(mat) > 0) return(as.vector(colSums(mat))) else
          return(numeric(0))
      }
      
      return(list(pgscore=do.call("mysum", lapply(results.list, function(x) as.vector(x$pgscore)))))
    }
  }
  #extract phenotypes
  phcovar <- parse.pheno.covar(pheno=NULL, covar=NULL, parsed=parsed, trace=2)
  pheno=phcovar$pheno
  pgscore=matrix(nrow=parsed$n, ncol=length(fit$lambda))
  
  if(is.null(parsed$extract)) {
    extract <- list(integer(0), integer(0))
    #col_skip_pos = integer(0); col_skip = integer(0)
  } else {
    #print(parsed$extract)
    extract <- selectregion(!parsed$extract)
    extract[[1]] <- extract[[1]] - 1
  }
  
  if(is.null(parsed$keep)) {
    keepbytes <- integer(0)
    keepoffset <- integer(0)
  } else {
    pos <- which(parsed$keep) - 1
    keepbytes <- floor(pos/4)
    keepoffset <- pos %% 4 * 2
  }
  #beta.id= which(extract.all) %in% which(parsed$extract) 
  #shrink=fit$shrink
  
  #r2=c()
  #AUC=c()
  library(pROC)
  for(i in 1:ncol(fit$beta)){
    #x1=x[,i,drop=F]
    #sd=fit$sd
    
    #pgs[,i]=(scale(test.data$G)/(sqrt(N-1)*sqrt(1-shrink))) %*% x1
    pgscore[,i]=pgs(fileName=fileName, N=parsed$N, P=parsed$P, shrink=fit$shrink, x=fit$beta[,i,drop=F], 
                    col_skip_pos=extract[[1]], col_skip=extract[[2]], keepbytes=keepbytes, keepoffset=keepoffset, 2)
    #r2[i]=cor(pgscore[,i], pheno, use="complete.obs") #remove missing values
    #AUC[i]=auc(roc(pheno, pgscore[,i]))
  }
  return(list(pgscore=pgscore
              #r2=r2,
              #AUC=AUC
  ))
}

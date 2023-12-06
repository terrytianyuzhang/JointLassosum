rm(list=ls()); gc()
options(stringsAsFactors = F)

#################
# code preps the data to be used in the JLS pipeline
# It is assumed that the input file are split by Chromosome
# On output, files with the essential information are generated for a common set of markers using the same effect/counted/A1 allele are generated.
# These files are set up to be used as the JLS input pipeline
#
# please read through the code to explore how the input.df should be set up
# all input files are required to be split by chromosome to facilitate using multiple cores for processing and reduce computing time. 
# The implemented version of multiple processing through doParallel does work on Windows operating systems.
# genotypes should be in the plink .bed, .bim, .fam format (standard plink1.9)
#
# plink1.9 is used as an external program and can be downloaded from https://www.cog-genomics.org/plink/1.9/

require(data.table)
require(doParallel)

plink="/data3/Software/Plink/plink"


###### Function
jlsDataPreparation <- function(chr,input.df,jls.dir,plink){
  # this function prepares the data to be be used in the JLS pipeline. 
  # datasets will be set up so that they have the same markers using the same ID and the same counted/A1/effect allele
  # in addition there will be a file with the original markers id mapping back to the consensus id
  # input:
  # chr - chromosome of interest
  # input.df - list in dataframe forward with all the essential information to create the unified JLS data
  # jls.dir - directory name to which the file should be written
  # plink - plink executable
  
  
  # set up variables to store data
  id.df=NULL
  DF.GWAS=list()  # to store the GWAS information
  
  # the next section to generate the consensus marker ID and it stores the GWAS information
  for(set in names(input.df)){
    # read the information on the SNP for each dataset
    if(tolower(input.df[[set]]$TYPE) == "gwas"){
      # this input data file is a GWA file and has a header
      tmp=fread(gsub("XX",chr,input.df[[set]]$inp.file),header=T,data.table=F)
      DF.GWAS[[set]]=tmp
    }else{
      # this is a plink file and we are reading the .bim file to get the SNP information
      tmp=fread(gsub("XX",chr,paste0(input.df[[set]]$inp.file,".bim")),header=F,data.table=F)
    }
    
    if(is.null(id.df)){
      # set up the consense id information based on the first dataset
      id.df=data.frame(ID=paste(gsub("chr","",toupper(tmp[,input.df[[set]]$CHR])),tmp[,input.df[[set]]$BP],tmp[,input.df[[set]]$CNTA],tmp[,input.df[[set]]$OTHA],sep=":"), # consensus ID
                       A1=tmp[,input.df[[set]]$CNTA])         # consensus counted/effect allele
      id.df[,paste0(set,".ID")]=tmp[,input.df[[set]]$ID]     # original ID 
    }else{
      # find the matches in the other datasets
      i12=match(id.df$ID,paste(gsub("CHR","",toupper(tmp[,input.df[[set]]$CHR])),tmp[,input.df[[set]]$BP],tmp[,input.df[[set]]$CNTA],tmp[,input.df[[set]]$OTHA],sep=":"))
      i21=match(id.df$ID,paste(gsub("CHR","",toupper(tmp[,input.df[[set]]$CHR])),tmp[,input.df[[set]]$BP],tmp[,input.df[[set]]$OTHA],tmp[,input.df[[set]]$CNTA],sep=":"))
      i12[is.na(i12)]=i21[is.na(i12)]
      id.df[,paste0(set,".ID")]=tmp[,input.df[[set]]$ID][i12]
      
    }
  }
  ### find the set of SNP respresented in all datasets
  sum(rowSums(is.na(id.df)) == 0)
  # select the data to be used, the markers that are available in all datasets
  id.df=id.df[rowSums(is.na(id.df)) == 0,]
  
  # remove duplicated entries
  id.df=id.df[!(id.df$ID %in% id.df$ID[duplicated(id.df$ID)]),]
  
  # save this to a cross-reference directory
  dir.create(paste0(jls.dir,"/XREF"), recursive = T,showWarnings = F)
  fwrite(id.df,paste0(jls.dir,"/XREF/marker-id-chr",chr,".xref"), row.names=F,col.names=T,quote=F,sep='\t')
  
  # in the GWAS data select the markers of interest
  for(set in names(DF.GWAS)){
    i.DF.GWAS=match(id.df[,paste0(set,".ID")],DF.GWAS[[set]][,input.df[[set]]$ID])
#    print(sum(is.na(i.DF.GWAS)))
    DF.GWAS[[set]]=DF.GWAS[[set]][i.DF.GWAS,]
  }
  
  ### process each of the datasets
  for(set in names(input.df)){
    if(tolower(input.df[[set]]$TYPE) == "gwas"){
      ### GWAS datasets
      df=data.frame(CHR=DF.GWAS[[set]][,input.df[[set]]$CHR], # chromosome
                    BP=DF.GWAS[[set]][,input.df[[set]]$BP],   # location
                    ID=id.df[,"ID"],                     # use consensus ID
                    A1=id.df[,"A1"])                      # use consensus counted/effect allele
      
      # effect size as odds ratio
      if(is.na(input.df[[set]]$B)){
        # odds ratio reported
        df[,"OR"]=DF.GWAS[[set]][,input.df[[set]]$OR]
      }else{
        # effect size reported, transform to odds ratio
        df[,"OR"]=exp(DF.GWAS[[set]][,input.df[[set]]$B])
      }
      # check if expressed on the consensus counted/effect allele
      # take 1/OR when not
      i.flip=which(df$A1 != DF.GWAS[[set]][,input.df[[set]]$CNTA])
      df[i.flip,"OR"]=1/df[i.flip,"OR"]
      
      # get the p-value
      df[,"P"]=DF.GWAS[[set]][,input.df[[set]]$P]
      
      # get the sample size
      if(is.numeric(input.df[[set]]$N)){
        # constant when reported as a number in the input.df
        df[,"N"]=input.df[[set]]$N
      }else{
        # use the value from the input file
        df[,"N"]=DF.GWAS[[set]][,input.df[[set]]$N]
      }
      
      # write the results to the outfile
      fwrite(df,gsub("XX",chr,input.df[[set]]$out.file),row.names=F,col.names=T,quote=F,sep="\t")
    }else{
      #### GENOTYPE DATASETS
      # write a file with the SNP to select to a temporary file      dir.create("Tmp",recursive = T,showWarnings = F)
      fwrite(id.df[,c(paste0(set,".ID"),"A1")],paste0("Tmp/",set,"-chr",chr,".a1"),row.names=F,col.names=F,quote=F,sep="\t")
      # write a file with samples to select to a temporary file
      fam=fread(paste0(gsub("XX",chr,input.df[[set]]$inp.file),".fam"),header=F,data.table=F)
      if("N" %in% colnames(input.df[[set]])){
        if(!is.na(input.df[[set]]$N)){
          fam=fam[1:min(nrow(fam),input.df[[set]]$N),]
        }
      }
      # the file with ids to keep from the plink data set
      fwrite(fam[,c("V1","V2")],paste0("Tmp/",set,"-chr",chr,".ind"),row.names=F,col.names=F,quote=F,sep="\t")
      ### create a plink genotype dataset with the selected markers and A1 as the counted allele
      plink.command=paste(plink,"--allow-no-sex","--nonfounders","--memory",10000,
                          "--bfile",gsub("XX",chr,input.df[[set]]$inp.file),
                          "--keep",paste0("Tmp/",set,"-chr",chr,".ind"),
                          "--extract",paste0("Tmp/",set,"-chr",chr,".a1"),
                          "--a1-allele",paste0("Tmp/",set,"-chr",chr,".a1"),
                          "--make-bed",
                          "--out",gsub("XX",chr,input.df[[set]]$out.file),
                          sep=" ")
      system(plink.command)
      # update the marker name to the consensus id
      bim=fread(paste0(gsub("XX",chr,input.df[[set]]$out.file),".bim"),header=F,data.table=F)
      i.id=match(bim$V2,id.df[,paste0(set,".ID")])
      bim$V2=id.df$ID[i.id]
      # write the bim file with the new marker information
      fwrite(bim,paste0(gsub("XX",chr,input.df[[set]]$out.file),".bim"),row.names=F,col.names=F,quote=F,sep="\t")
      
      # clean the Tmp directory
      system(paste("rm",paste0("Tmp/",set,"-chr",chr,".a1"),paste0("Tmp/",set,"-chr",chr,".ind")))
    }
  }
  
  # get back to main
  return(data.frame(chr,n.marker=nrow(id.df)))
}

################# END FUNCTION

# specify the master location for the jls input data files
jls.dir="JLS"
dir.create(jls.dir,recursive = T, showWarnings = F)
# create a temporary directory
dir.create("Tmp",recursive = T, showWarnings = F)


input.df=list()
# Data for which PRS needs to be determine, this data is in .bed, .bim, .fam format and organized by chromosome, these do not have header so V1-V6 is used
# REQUIRED
input.df[["PRS"]]=data.frame(TYPE="genotype",              # type of dataset (genotype or gwas)
                             CHR="V1",                     # column name in input with chromosome number
                             BP="V4",                      # column name in input with location
                             ID="V2",                     # column name in input with marker id
                             CNTA="V5",                    # column name in input with effect or counted allele
                             OTHA="V6",                    # column name in input with effect or counted allele
                             inp.file="PRS/byChr/PRS-hg38-chrXX", # input location and file name, using XX for the chromosome (do not use a file extension)
                             out.file=paste0(jls.dir,"/PRS/PRS-hg38-chrXX")) # output location and file name, using XX for the chromosome (do not use a file extension)

# first GWAS
# REQUIRED
input.df[["GWAS1"]]=data.frame(TYPE="gwas",              # type of dataset (genotype or gwas)
                               CHR="Chromosome",            # column name in input with chromosome number
                               BP="Position",               # column name in input with location
                               ID="MarkerName",            # column name in input with marker id
                               CNTA="Effect_allele",         # column name in input with effect or counted allele
                               OTHA="Non_Effect_allele",    # column name in input with effect or counted allele
                               B = "Beta",                  # column name in input with effect size (use NA when not available, either B or OR must be specified)
                               OR = NA,                     # column name in input with odd ratio (use NA when not available, either B or OR must be specified)
                               P = "Pvalue",                # column name in input with p-value
                               N = 10000,                   # column name in input with number of samples or the number of samples to be used for all markers
                               inp.file="GWAS/EUR/byChr/EUR-hg38-chrXX.gwas.gz", # input location and file name, using XX for the chromosome number
                               out.file=paste0(jls.dir,"/GWAS1/EUR-hg38-chrXX.gwas.gz")) # output location and file name, using XX for the chromosome
                               
# second GWAS
# REQUIRED
input.df[["GWAS2"]]=data.frame(TYPE="gwas",              # type of dataset (genotype or gwas)
                               CHR="CHR",                    # column name in input with chromosome number
                               BP="BP",                      # column name in input with location
                               ID="SNP",                     # column name in input with marker id
                               CNTA="A1",                    # column name in input with effect or counted allele
                               OTHA="A2",                    # column name in input with effect or counted allele
                               B = NA,                       # column name in input with effect size (use NA when not available, either B or OR must be specified)
                               OR = "OR",                    # column name in input with odd ratio (use NA when not available, either B or OR must be specified)
                               P = "P",                      # column name in input with p-value
                               N = "Neff",                   # column name in input with number of samples or the number of samples to be used for all markers
                               inp.file="GWAS/EAS/byChr/EAS-hg38-chrXX.gwas.gz", # input location and file name, using XX for the chromosome number
                               out.file=paste0(jls.dir,"/GWAS2/EAS-hg38-chrXX.gwas.gz")) # output location and file name, using XX for the chromosome

# Genotype data for the first REFERENCE population matching in ancestry to the first GWAS,
# it can be either genotypes from random set of samples from the same ancestry, or one can use the synthetic populations as the set
# this is not used when you use the synthetic population for the parameter tuninig. 
# Often these reference populations are quite large and not all samples are needed. By specifying N to be less than the total number of samples
# computations can be speeded up.
# the files are assumed to be split by chromosome and are in the plink .bed, .bim, .fam format
# REQUIRED
input.df[["REF1"]]=data.frame(TYPE="genotype",              # type of dataset (genotype or gwas)
                              CHR="V1",                    # column name in input with chromosome number
                              BP="V4",                      # column name in input with location
                              ID="V2",                    # column name in input with marker id
                              CNTA="V5",                    # column name in input with effect or counted allele
                              OTHA="V6",                    # column name in input with effect or counted allele
                              N= 4000,                     # specify the number of samples from this population to use (see comment above)
                              inp.file="../Data/EUR-Synthetic/byChr/EUR-chrXX", # input location and file name, using XX for the chromosome (do not use a file extension)
                              out.file=paste0(jls.dir,"/REF1/EUR-hg38-chrXX")) # output location and file name, using XX for the chromosome (do not use a file extension)

# Genotype data for the second REFERENCE population matching in ancestry to the second GWAS,
# it can be either genotypes from a random set of samples from the same ancestry, or one can use the synthetic populations as the set
# this is not used when you use the synthetic population for the parameter tuninig. 
# Often these reference populations are quite large and not all samples are needed. By specifying N to be less than the total number of samples
# computations can be speeded up.
# the files are assumed to be split by chromosome and are in the plink .bed, .bim, .fam format
# REQUIRED
input.df[["REF2"]]=data.frame(TYPE="genotype",              # type of dataset (genotype or gwas)
                              CHR="V1",                    # column name in input with chromosome number
                              BP="V4",                      # column name in input with location
                              ID="V2",                    # column name in input with marker id
                              CNTA="V5",                    # column name in input with effect or counted allele
                              OTHA="V6",                    # column name in input with effect or counted allele
                              N= 4000,                     # specify the number of samples from this population to use (see comment above)
                              inp.file="../Data/EAS-Synthetic/byChr/EAS-chrXX", # input location and file name, using XX for the chromosome (do not use a file extension)
                              out.file=paste0(jls.dir,"/REF2/EAS-hg38-chrXX")) # output location and file name, using XX for the chromosome (do not use a file extension)

# Genotype data for the first SYNTHETIC population matching in ancestry to the first GWAS,
# the files are assumed to be split by chromosome and are in the plink .bed, .bim, .fam format
# when not used set input.df[["SYN1"]]=NULL
# REQUIRED only when using SYNTHETIC for parameter tuning
input.df[["SYN1"]]=data.frame(TYPE="genotype",              # type of dataset (genotype or gwas)
                              CHR="V1",                    # column name in input with chromosome number
                              BP="V4",                      # column name in input with location
                              ID="V2",                    # column name in input with marker id
                              CNTA="V5",                    # column name in input with effect or counted allele
                              OTHA="V6",                    # column name in input with effect or counted allele
                              inp.file="../Data/EUR-Synthetic/byChr/EUR-chrXX", # input location and file name, using XX for the chromosome (do not use a file extension)
                              out.file=paste0(jls.dir,"/SYN1/EUR-hg38-chrXX")) # output location and file name, using XX for the chromosome (do not use a file extension)

# Genotype data for the second SYNTHETIC population matching in ancestry to the first GWAS,
# the files are assumed to be split by chromosome and are in the plink .bed, .bim, .fam format
# when not used set input.df[["SYN1"]]=NULL
# REQUIRED only when using SYNTHETIC for parameter tuning
input.df[["SYN2"]]=data.frame(TYPE="genotype",              # type of dataset (genotype or gwas)
                              CHR="V1",                    # column name in input with chromosome number
                              BP="V4",                      # column name in input with location
                              ID="V2",                    # column name in input with marker id
                              CNTA="V5",                    # column name in input with effect or counted allele
                              OTHA="V6",                    # column name in input with effect or counted allele
                              inp.file="../Data/EAS-Synthetic/byChr/EAS-chrXX", # input location and file name, using XX for the chromosome (do not use a file extension)
                              out.file=paste0(jls.dir,"/SYN2/EAS-hg38-chrXX")) # output location and file name, using XX for the chromosome (do not use a file extension)


# create all the required output directories
for(set in names(input.df)){
  dir=unlist(strsplit(gsub(paste0(jls.dir,"/"),"",input.df[[set]]$out.file),"/")); dir=dir[-length(dir)]
  dir.create(paste(c(jls.dir,dir),collapse="/"),recursive = T,showWarnings = F)
}

#### prepare the data for JLS by chr for all the pieces of information
system.time(re<-unlist(mclapply(1:22,jlsDataPreparation,input.df=input.df,jls.dir=jls.dir,plink=plink,mc.cores=22,mc.preschedule = F, mc.silent=T)))

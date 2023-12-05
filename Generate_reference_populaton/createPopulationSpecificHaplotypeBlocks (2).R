rm(list=ls()); gc()
options(stringsAsFactors = F)

#####
### create haplotypes for haplotypes from the selected 1000G genotypes from step 1.
### the user needs to specify the super populations to use (EUR, AFR, ASI (for EAS)) (easiest is to use the super population and population information from the previous step)
### build is assumed to be GRCh38/hg38
### selected markers are represented in haplotype block from each of the super-populations

### sofware used
### set link to external software
plink2="/data3/Software/Plink2/plink2"
bcftools="/data3/Software/SAMTOOLS/bcftools-1.16/bcftools"

# set directories
tmp.dir="Tmp/"                                                  # used for intermediate data files
data.dir="Data/"                                                # used for the final data files

# libraries used
require(data.table)
require(doParallel)
require(lassosum)            # contains haplotype block information from Berisa and Pickrell (2015, PMID: 26395773)

###### FUNCTION USED
# function used in the code
createHaplotypes <- function(i.block,pop,blocks,haplo.dir,pop.haplo.dir){
  ### this function creates the haplotype file for each region
  # i.block - block of interest
  # pop - super population
  # blocks - ld-blocks for the super-population
  # haplo.dir - directory in which the combined haplotype files from step (1) are located
  # pop.haplo.dir - directory in which the population and ld block specific haplotypes are stored
  region=paste0(blocks$chr[i.block],":",blocks$start[i.block],"-",blocks$stop[i.block]) # create the region
  
  # select the combined haplotypes in the region
  bcf.command=paste(bcftools,"view",
                    paste0(haplo.dir,"/EUR-AFR-EAS-",blocks$chr[i.block],".vcf.gz"),
                    paste("--samples-file",paste0(pop.haplo.dir,"/",pop,".samples")), # samples
                    paste("--regions",region),                                                                         # region
                    paste("--output",paste0("Tmp/",pop,"-",region,".vcf.gz")))                                         # output to a temporary file
  system(bcf.command)
  
  # read the combined haplotypes and split them
  n=length(grep('##',readLines(paste0("Tmp/",pop,"-",region,".vcf.gz"),n = 150)))
  gnt=fread(paste0("Tmp/",pop,"-",region,".vcf.gz"),header=T,data.table=F,skip = n,sep="\t")
  rownames(gnt)=gnt$ID
  map=gnt[,c("#CHROM","ID","POS","REF","ALT")]
  gnt=as.matrix(gnt[,-c(1:9)])            # remove the information part
  gnt=gsub("|"," ",gnt,fixed = T)         # change the delimiter in the haplotype calls for a sample to a space
  fwrite(as.data.frame(gnt),paste0("Tmp/",pop,"-",region,".tmp.gz"),row.names=F,col.names=F,quote=F,sep=" ")  # write this information to a file so that it can be read back as space delimited, 
  
  # read the split haplotypes and process them                                                                                                            # splitting the combined haplotype
  haplos=fread(paste0("Tmp/",pop,"-",region,".tmp.gz"),header=F,data.table=F,sep=" ")
  colnames(haplos)=paste0("H",1:ncol(haplos))                                              # give names to the columns
  rownames(haplos)=map$ID                                                                  # give rownames to the haplotype blocks
  # write the haplotype block
  fwrite(map,paste0(pop.haplo.dir,"/",pop,"-",region,".map"),row.names=F,col.names=T,quote=F,sep="\t")
  fwrite(haplos,paste0(pop.haplo.dir,"/",pop,"-",region,".haplo.gz"),row.names=T,col.names=T,quote=F,sep="\t")
  
  # perform clean-up
  system(paste("rm",paste0("Tmp/",pop,"-",region,".vcf.gz"),paste0("Tmp/",pop,"-",region,".tmp.gz")))
  
  return(data.frame(i.block,region,n.marker=nrow(gnt)))
}
#### END FUNCTION

#### load the super-population information from step (1)
load("Data/super.populations.RData")
# get the population information
pop.BP=c("EUR","AFR","ASN"); names(pop.BP)=c("EUR","AFR","EAS")  #  #  Berisa and Pickrell used ASN instead of EAS to designate the East Asians


##### EUROPEAN SUPER-POPULATION HAPLOTYPE FILES
pop="EUR"
haplo.dir=paste0("Data/","1000G-Selected-Haplotypes")               # directory with the combined haplotypes for all populations
pop.haplo.dir=paste0("Data/","1000G-Selected-HaplotypeBlocks-",pop) # population specific haplotype directory, will be created
# create directories
dir.create(pop.haplo.dir,showWarnings = F, recursive = T)

# write a list of samples for the population of interest
fwrite(as.data.frame(super.pop[[pop]]$samples),paste0(pop.haplo.dir,"/",pop,".samples"),row.names=F,col.names=F,quote=F,sep="\t")

# read the boundaries of the haplotype block from lassosum
blocks=fread(system.file(paste0("data/Berisa.",pop.BP[pop],".hg38", ".bed"),package="lassosum"), header=T,data.table=F)
# since there is overlap reduce the stop location by one base pair
blocks$stop=blocks$stop-1
# create the haplotype block files
system.time(re<-rbindlist(mclapply(sample(nrow(blocks)),createHaplotypes,pop=pop,blocks=blocks,
                                   haplo.dir=haplo.dir,pop.haplo.dir=pop.haplo.dir,mc.cores=32,mc.preschedule = F)))
sum(re$n.marker)
re=re[order(re$i.block),]
fwrite(re,paste0(pop.haplo.dir,"/",pop,"-haplo-block-information.txt"),row.names=F,col.names=T,quote=F,sep="\t")

##### AFRICAN SUPER-POPULATION HAPLOTYPE FILES
pop="AFR"
haplo.dir=paste0("Data/","1000G-Selected-Haplotypes")               # directory with the combined haplotypes for all populations
pop.haplo.dir=paste0("Data/","1000G-Selected-HaplotypeBlocks-",pop) # population specific haplotype directory, will be created
# create directories
dir.create(pop.haplo.dir,showWarnings = F, recursive = T)

# write a list of samples for the population of interest
fwrite(as.data.frame(super.pop[[pop]]$samples),paste0(pop.haplo.dir,"/",pop,".samples"),row.names=F,col.names=F,quote=F,sep="\t")

# read the boundaries of the haplotype block from lassosum
blocks=fread(system.file(paste0("data/Berisa.",pop.BP[pop],".hg38", ".bed"),package="lassosum"), header=T,data.table=F)
# since there is overlap reduce the stop location by one base pair
blocks$stop=blocks$stop-1
# create the haplotype block files
system.time(re<-rbindlist(mclapply(sample(nrow(blocks)),createHaplotypes,pop=pop,blocks=blocks,
                                   haplo.dir=haplo.dir,pop.haplo.dir=pop.haplo.dir,mc.cores=32,mc.preschedule = F)))
sum(re$n.marker)
re=re[order(re$i.block),]
fwrite(re,paste0(pop.haplo.dir,"/",pop,"-haplo-block-information.txt"),row.names=F,col.names=T,quote=F,sep="\t")

##### EAST ASIAN SUPER-POPULATION HAPLOTYPE FILES
pop="EAS"
haplo.dir=paste0("Data/","1000G-Selected-Haplotypes")               # directory with the combined haplotypes for all populations
pop.haplo.dir=paste0("Data/","1000G-Selected-HaplotypeBlocks-",pop) # population specific haplotype directory, will be created
# create directories
dir.create(pop.haplo.dir,showWarnings = F, recursive = T)

# write a list of samples for the population of interest
fwrite(as.data.frame(super.pop[[pop]]$samples),paste0(pop.haplo.dir,"/",pop,".samples"),row.names=F,col.names=F,quote=F,sep="\t")

# read the boundaries of the haplotype block from lassosum
blocks=fread(system.file(paste0("data/Berisa.",pop.BP[pop],".hg38", ".bed"),package="lassosum"), header=T,data.table=F)
# since there is overlap reduce the stop location by one base pair
blocks$stop=blocks$stop-1
# create the haplotype block files
system.time(re<-rbindlist(mclapply(sample(nrow(blocks)),createHaplotypes,pop=pop,blocks=blocks,
                                   haplo.dir=haplo.dir,pop.haplo.dir=pop.haplo.dir,mc.cores=32,mc.preschedule = F)))
sum(re$n.marker)
re=re[order(re$i.block),]
fwrite(re,paste0(pop.haplo.dir,"/",pop,"-haplo-block-information.txt"),row.names=F,col.names=T,quote=F,sep="\t")

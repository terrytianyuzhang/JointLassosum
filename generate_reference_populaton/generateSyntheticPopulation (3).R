rm(list=ls()); gc()
options(stringsAsFactors = F)

################
# this code generates the synthetic populations based on the haplotypes that were created in step (2)
# data is processed by blocks and can easily be run in parallel
# pop - super population
# N - number of samples to generate for the synthethic population

require(data.table)
require(doParallel)
require(snpStats)

plink="/data3/Software/Plink/plink"
plink2="/data3/Software/Plink2/plink2"

##### functions used
synthethicGenotypesByBlock <- function(region,pop,N,haplo.dir,by.block.dir){
  # this function generates the synthetic genotypes by block
  # region - region to generate
  # pop - population
  # N - number of samples to generate
  # haplo.dir - directory where the haplotypes are stored
  # by,block.dir - directory for the by block genotypes
  
  # read the map and haplotype information
  map=fread(paste0(haplo.dir,"/",pop,"-",region,".map"),header=T,data.table=F)
  haplos=fread(paste0(haplo.dir,"/",pop,"-",region,".haplo.gz"),header=T,data.table=F)
  rownames(haplos)=haplos$V1               # use marker names in V1 as the rownames
  haplos=t(as.matrix(haplos[,-1]))         # remove the marker name column and transpose the haplotype matrix for easy processing
  
  # randomly pick two haplotypes for each sample to generate
  haplo1=sample(rownames(haplos),N,replace=T)         # first haplotype
  haplo2=sample(rownames(haplos),N,replace=T)         # second haplotype
  while(sum(haplo1 == haplo2)){                       # don't allow haplo1 to be the same as haplo2
    i=which(haplo1 == haplo2)
    haplo2[i]=sample(rownames(haplos),length(i),replace=T) # replace the ones that are the same with a new one for haplo2
  }
  
  # generate the genotypes by adding the allele count from each haplotype
  gnt=haplos[haplo1,]+haplos[haplo2,]        
  rownames(gnt)=paste0(pop,1:N) # generate sample names
  
  # write the genotypes to a plink binary format file using write.plink in snpStats
  write.plink(paste0(by.block.dir,"/",pop,"-",region),
              snp.major=T,
              snps=as(gnt,"SnpMatrix"),
              pedigree=rownames(gnt),
              id=rownames(gnt),
              phenotype=-9,
              chromosome=gsub("chr","",map$`#CHROM`),
              genetic.distance=rep(0,ncol(gnt)),
              position=map$POS,
              allele.1=map$ALT,
              allele.2=map$REF,
              na.code=0)
  # write a file with the counted allele
  fwrite(map[,c("ID","ALT")],paste0(by.block.dir,"/",pop,"-",region,".A1"),row.names=F,col.names=F,quote=F,sep="\t")
  
  return(region)
  
}


mergeBlocksIntoChromosomes <- function(chr,pop,blocks,plink,by.block.dir,by.chr.dir){
  # merge the by block genotypes to by chromosome
  # chr - chromosome
  # pop - population
  # blocks -ld blocking information
  # by.block.dir - location of the genotypes by blocks
  # by.chr.dir - location of the genotypes by chromsome
  
  sel.blocks=blocks[grep(paste0("chr",chr,":"),blocks$region),"region"]
  sel.blocks=paste0(by.block.dir,"/",pop,"-",sel.blocks)  # turn them into the location of the plink genotype files
  # create directory to store the genotypes by chromosome
  # store the selected blocks to be combined
  fwrite(as.data.frame(sel.blocks[-1]),paste0(by.chr.dir,"/",pop,"-blocks-chr",chr,".merge"),row.names=F,col.names=F,quote=F,sep="\t")
  # merge the blocks into chromosomes
  plink.command=paste(plink,"--allow-no-sex","--nonfounders","--memory",20000,
                      "--bfile",sel.blocks[1],
                      "--keep-allele-order",
                      "--merge-list",paste0(by.chr.dir,"/",pop,"-blocks-chr",chr,".merge"),
                      "--make-bed",
                      "--out",paste0(by.chr.dir,"/",pop,"-chr",chr),
                      sep=" ")
  system(plink.command)
  # collect information on A1 allele
  a1.allele=NULL
  for(file in sel.blocks){
    a1.allele=rbind(a1.allele,fread(paste0(file,".A1"),header=F,data.table=F))
  } 
  fwrite(a1.allele,paste0(by.chr.dir,"/",pop,"-chr",chr,".A1"),row.names=F,col.names=F,quote=F,sep="\t")
  
  system(paste("rm",paste0(by.chr.dir,"/",pop,"-chr",chr,".nosex")))
  return(data.frame(chr,n.block=length(sel.blocks),n.marker=nrow(a1.allele)))
  
}

# generate synthetics for each population of interest
pop="EUR"    ### population
N=20000       ### number of random samples to create

haplo.dir=paste0("Data/1000G-Selected-HaplotypeBlocks-",pop)  # population specific haplotype directory
by.block.dir=paste0("Data/",pop,"-Synthetic/byBlock")         # population specific genotype by block directory
by.chr.dir=paste0("Data/",pop,"-Synthetic/byChr")             # population specific genotype by chromosome directory
# create the two genotype directories
dir.create(by.block.dir,recursive = T, showWarnings = F)         
dir.create(by.chr.dir,recursive = T, showWarnings = F)           

### read the information on the haplotype blocks
blocks=fread(paste0(haplo.dir,"/",pop,"-haplo-block-information.txt"),header=T,data.table=F)
# generate genotypes for each block
system.time(re<-unlist(mclapply(blocks$region,synthethicGenotypesByBlock,pop=pop,N=N,
                                haplo.dir=haplo.dir,by.block.dir=by.block.dir,mc.cores=64,mc.preschedule = F, mc.silent=T)))

# merge the block genotypes into chromsomes
system.time(re<-rbindlist(mclapply(1:22,mergeBlocksIntoChromosomes,pop = pop,blocks = blocks,
                                   plink = plink,by.block.dir = by.block.dir,by.chr.dir = by.chr.dir,
                                   mc.cores = 22,mc.preschedule = T,mc.silent = T)))
# this will give warning message about same position markers, this can be ignored


pop="AFR"    ### population
N=4000       ### number of random samples to create

haplo.dir=paste0("Data/1000G-Selected-HaplotypeBlocks-",pop)  # population specific haplotype directory
by.block.dir=paste0("Data/",pop,"-Synthetic/byBlock")         # population specific genotype by block directory
by.chr.dir=paste0("Data/",pop,"-Synthetic/byChr")             # population specific genotype by chromosome directory
# create the two genotype directories
dir.create(by.block.dir,recursive = T, showWarnings = F)         
dir.create(by.chr.dir,recursive = T, showWarnings = F)           

### read the information on the haplotype blocks
blocks=fread(paste0(haplo.dir,"/",pop,"-haplo-block-information.txt"),header=T,data.table=F)
# generate genotypes for each block
system.time(re<-unlist(mclapply(blocks$region,synthethicGenotypesByBlock,pop=pop,N=N,
                                haplo.dir=haplo.dir,by.block.dir=by.block.dir,mc.cores=64,mc.preschedule = F, mc.silent=T)))

# merge the block genotypes into chromsomes
system.time(re<-rbindlist(mclapply(1:22,mergeBlocksIntoChromosomes,pop = pop,blocks = blocks,
                                   plink = plink,by.block.dir = by.block.dir,by.chr.dir = by.chr.dir,
                                   mc.cores = 22,mc.preschedule = T,mc.silent = T)))
# this will give warning message about same position markers, this can be ignored


pop="EAS"    ### population
N=10000      ### number of random samples to create

haplo.dir=paste0("Data/1000G-Selected-HaplotypeBlocks-",pop)  # population specific haplotype directory
by.block.dir=paste0("Data/",pop,"-Synthetic/byBlock")         # population specific genotype by block directory
by.chr.dir=paste0("Data/",pop,"-Synthetic/byChr")             # population specific genotype by chromosome directory
# create the two genotype directories
dir.create(by.block.dir,recursive = T, showWarnings = F)         
dir.create(by.chr.dir,recursive = T, showWarnings = F)           

### read the information on the haplotype blocks
blocks=fread(paste0(haplo.dir,"/",pop,"-haplo-block-information.txt"),header=T,data.table=F)
# generate genotypes for each block
system.time(re<-unlist(mclapply(blocks$region,synthethicGenotypesByBlock,pop=pop,N=N,
                                haplo.dir=haplo.dir,by.block.dir=by.block.dir,mc.cores=64,mc.preschedule = F, mc.silent=T)))

# merge the block genotypes into chromsomes
system.time(re<-rbindlist(mclapply(1:22,mergeBlocksIntoChromosomes,pop = pop,blocks = blocks,
                                   plink = plink,by.block.dir = by.block.dir,by.chr.dir = by.chr.dir,
                                   mc.cores = 22,mc.preschedule = T,mc.silent = T)))
# this will give warning message about same position markers, this can be ignored

rm(list=ls()); gc()
options(stringsAsFactors = F)

### selectDataFrom1000G
# this code selects variants and genotypes from the 1000G haplotypes build hg38
# software needed:  - bcftools (http://www.htslib.org/download/)
#                   - plink2 (https://www.cog-genomics.org/plink/2.0/#recent)
#
# data needed: - copy of the phased haplotypes from 1000G (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/)
#              - sample information file (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/20130606_g1k_3202_samples_ped_population.txt)   
#
# at the end of this code you will have haplotypes for selected variants and selected samples.
#
# JLS only uses the autosomes. The haplotype blocks that it relies on and are defined in Berisa and Pickrell (2015, PMID: 26395773) are only published for the autosomes.

require(data.table)
require(doParallel)

### set link to external software
plink2="/data3/Software/Plink2/plink2"
bcftools="/data3/Software/SAMTOOLS/bcftools-1.16/bcftools"

# set directories
samples.dir="/data3/DownLoadedData/1000G/20201028_3202_phased/" # set the location of the sample file
orig.haplo.dir="/data3/DownLoadedData/1000G/20201028_3202_phased/"   # set the location of the downloaded haplotype files
tmp.dir="Tmp/"                                                  # used for intermediate data files
data.dir="Data/"                                                # used for the final data files

# create directories that are being used
sel.haplo.dir=paste0(data.dir,"1000G-Selected-Haplotypes")
dir.create(sel.haplo.dir,recursive = T, showWarnings = F)


# specify the populations to use (can specify multipe populations for each super population
# supported super population are EUR, AFR, EAS))
super.pop=list()
super.pop[["EUR"]]=list()
super.pop[["EUR"]]$pop=c("CEU")
super.pop[["AFR"]]=list()
super.pop[["AFR"]]$pop=c("YRI")
super.pop[["EAS"]]=list()
super.pop[["EAS"]]$pop=c("CHB","CHS")

# find the samples for each super.pop
samples=fread("/data3/DownLoadedData/1000G/Phase3-Samples/20130606_g1k_3202_samples_ped_population.txt",header=T,data.table=F)
select.ids=NULL
for(pop in names(super.pop)){
  super.pop[[pop]]$samples=samples$SampleID[samples$Population %in% super.pop[[pop]]$pop]
  select.ids=c(select.ids,super.pop[[pop]]$samples)
}
# write a list of sample names
select.id.file=paste0(sel.haplo.dir,"/selected.ids")
fwrite(as.data.frame(select.ids),select.id.file,row.names = F,col.names=F,quote=F,sep="\t")
# save the super populations for future steps
save(super.pop,file="Data/super.populations.RData")

########################################################
########################################################
########################################################
# SETTING UP THE FILTERS TO USE FOR SELECTING SAMPLES AND VARIANTS FROM THE 1000G FILES
########################################################
########################################################
########################################################
# select the data to use from the 1000G haplotype files
min.maf=0.01     # minimum minor allele frequency in each of the super populations
min.phwe=0.005
bcf.options="-m2 -M2 -v snps"     # bcftools command to filter for biallelic SNP, see bcftools documentation for other options
                                  # http://samtools.github.io/bcftools/bcftools.html#view
bcf.select=paste("--samples-file",select.id.file)  # define the samples to extract from the vcf file
# set up the filters, this is to include variants
maf.filter=paste(paste(paste0("INFO/AF_",names(super.pop),"_unrel"),">=",min.maf,"&",paste0("INFO/AF_",names(super.pop),"_unrel"),"<=",1-min.maf),collapse=" & ")
hwe.filter=paste(paste(paste0("INFO/HWE_",names(super.pop)),">=",min.phwe),collapse=" & ")

# combine the options used for view
view.options=paste(bcf.options,bcf.select,paste0("--include ","'",paste(maf.filter,hwe.filter,sep=" & "),"'"))

### name output file do not suppy a file extension, this will be placed in directory sel.haplo.dir, one for each chromosome, appended with -chrXX.vcf.gz
outfile=paste0(sel.haplo.dir,"/",paste(names(super.pop),collapse="-"))


select1000GHaplotypes <- function(chr,bcftools,view.options,orig.haplo.dir,outfile){
  # this function selects the haplotypes for the super populations from the 1000G haplotype files 
  # CCDG_14151_B01_GRM_WGS_2020-08-05_chrXX.filtered.shapeit2-duohmm-phased.vcf.gz, where XX is the chromosome names (1:22)
  # 
  # chr: chromosome number
  # bcftools: bcftools executable
  # view.options: string of options to use when running bcftools with command view
  # haplo.dir: location of the downloaded haplotype files
  # outfile: name of the file to which the selected variants and samples are written
  # it can be run in parallel one for each chromosome
  
  bcf.command=paste(bcftools,"view",
                    paste0(orig.haplo.dir,"CCDG_14151_B01_GRM_WGS_2020-08-05_chr",chr,".filtered.shapeit2-duohmm-phased.vcf.gz"),
                    view.options,
                    "-o",paste0(outfile,"-chr",chr,".vcf.gz"))
  system(bcf.command)
  # tabix the file
  system(paste("tabix -f -p vcf",paste0(outfile,"-chr",chr,".vcf.gz")))
  
  return(chr)
}

### run for all chromosomes
mclapply(1:22,select1000GHaplotypes,bcftools = bcftools,view.options = view.options,orig.haplo.dir = orig.haplo.dir,outfile = outfile,mc.cores=24,mc.preschedule = F)


##### NOT NEEDED
### pre-tabix the downloaded files to make sure that they are up to date
#runTabix <- function(chr,orig.haplo.dir){
#  system(paste("tabix -f -p vcf",paste0(orig.haplo.dir,"CCDG_14151_B01_GRM_WGS_2020-08-05_chr",chr,".filtered.shapeit2-duohmm-phased.vcf.gz")))
#  return(chr)
#}
#mclapply(1:22,runTabix,orig.haplo.dir = orig.haplo.dir,mc.cores=24,mc.preschedule=T)

### arguments
  args = commandArgs(trailingOnly=TRUE)
  fai.fn <- args[1]
  focalFile.fn <- args[2]
  species <- args[3]
  nJobs <- as.numeric(args[4])

### libraries
  library(data.table)
  library(foreach)

### load data
  #fai.fn="/Users/alanbergland/Documents/GitHub/DESTv3/snpCalling_dev/scatter_gather_annotate/holo_dmel_6.12.fa.fai"
  #focalFile.fn="/Users/alanbergland/Documents/GitHub/DESTv3/examples/mapping/focalFile"
  #species <- "snail"
  #nJobs <- 5000

  fai.fn="/scratch/aob2x/snail/refGenome/GCF_037325665.1_US_GU_Lsax_2.0_genomic.fa.fai"
  focalFile.fn="/scratch/aob2x/snail/focalFile"
  nJobs=2000
  species="snail"
  print(fai.fn)
  print(focalFile.fn)

  fai <- fread(fai.fn, header=F)
  focalFile <- fread(focalFile.fn, header=F)

### subset to species
  setkey(fai, V1)
  use <- fai[J(focalFile[V1==species]$V2)]
  use <- use[!is.na(V2)]
  use[,cLen:=cumsum(V2)]
  use[,pLen:=V2/max(cLen)]
  use[,nJobs_chr:=floor(nJobs*pLen)]
  sum(use$nJobs_chr)

### iterate across chromosomes
  wins <- foreach(i=c(1:dim(use)[1]), .combine="rbind")%do%{
    # i <- 1
    data.table(chr=use$V1[i],
               start=floor(seq(from=1, to=use$V2[i], length.out=(1+use$nJobs_chr[i])))[-(1+use$nJobs_chr[i])],
               stop= floor(seq(from=1, to=use$V2[i], length.out=(1+use$nJobs_chr[i])))[-1]-1)
  }
  dim(wins)

### export
  write.table(wins, file=paste(fai.fn, species, nJobs, "jobs", sep="."), sep=",", quote=F, row.names=F, col.names=F)

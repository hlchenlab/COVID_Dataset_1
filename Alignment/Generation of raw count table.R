#################################

##Author Name: Conor Cremin
##Title: Parsing of STAR Output into dataframe format
##Dataset:COVOD-19

#################################

# Step 1: Set parameters:

gseid= "GSE156006"
files = matrix(unlist(strsplit(as.character(list.files("/path/to/directory for STAR output",pattern = "ReadsPerGene.out.tab")), "Reads")),ncol=2, byro=T)[,1]
genecounts = "ReadsPerGene.out.tab"
logname = "Log.final.out"
filedir = paste("/path/to/directory for STAR output","", sep="/") # Specify file path

# Step 2: Counting Loop formula:

Ns = list()
i = 1
for( n in files ){
  N = list()
  countfile = paste0(filedir,n, genecounts)
  logfile = paste0(filedir,n, logname)
  
  if( file.exists(countfile) ) {
    print(countfile)
    counts =  read.table(countfile)
    
    log1 =read.table(logfile, sep="\t", nrows=6)
    log2 =read.table(logfile, sep="\t", skip=8, nrows=14)
    log3 =read.table(logfile, sep="\t", skip=23, nrows=4)
    log4 =read.table(logfile, sep="\t", skip=28, nrows=3)
    
    N$mapinfo = rbind(log1,log2,log3,log4)
    N$unmapped =  counts[1,]
    N$multimapping = counts[2,]
    N$noFeature =   counts[3,]
    N$ambiguous = counts[4,]
    N$length = dim(counts)[1]-4
    N$genes = counts[ (1:N$length)+4,1]
    N$counts1 = counts[ (1:N$length)+4,2]  #Unstranded
    N$counts2 = counts[ (1:N$length)+4,3]  #Forward
    N$counts3 = counts[ (1:N$length)+4,4]  #Reverse
    
  } else {
    N$counts3 = rep(0, length(files) ) ## Change if appropriate counts are in col 3/4
  }
  if( i > 1  ){
    counts_exp = cbind(counts_exp, N$counts3) ## Change if appropriate counts are in col 3/4
  } else {
    counts_exp = N$counts3  ## Change if appropriate counts are in col 3/4
  }
  Ns[[i]] = N
  print(i)
  i = i + 1
}

# Step 3: Format Row Names to remove the decimal from ENSG names and Save:

genes = list()
colnames(counts_exp) = as.character(files) 
if( !require("stringr")){
  BiocManager::install("stringr")
}
for (i in 1:length(N$genes)) {
  if(str_detect(N$genes[i],".")){
    genes[i] = matrix(unlist(strsplit( as.character(N$genes[i]), "\\.") ) , ncol=2, byro=T)[,1]
  } else {
    genes[i] = as.character(N$genes[i])
  }
}
genes = do.call(rbind.data.frame,genes)
rownames(counts_exp) = as.character(genes[,1])

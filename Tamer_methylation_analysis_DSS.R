library(DSS)


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) #take commandline args
if(length(args) < 5) {
  args <- c("--help")
}

S1_R1 <- as.character(unlist(args[1])) #sample 1 file
S1_R2 <- as.character(unlist(args[2])) #sample 1 file

S2_R1 <- as.character(unlist(args[3])) #sample 2 file
S2_R2 <- as.character(unlist(args[4])) #sample 2 file

output_rda_file <- as.character(unlist(args[5]))
output_DML_file <- as.character(unlist(args[6]))
output_DMR_file <- as.character(unlist(args[7]))

require(bssseq)
tab5rows <- read.table(sample_1, header = F,sep="\t", nrows = 20)
classes <- sapply(tab5rows, class)
Sample_1 <- read.table(sample_1, header = F,sep="\t", colClasses = classes)
Sample_2 <- read.table(sample_2, header = F,sep="\t", colClasses = classes)
colnames(Sample_1) <- c("chr","pos","N","X")
colnames(Sample_2) <- c("chr","pos","N","X")
sub_Sample_1 <- subset(Sample_1,N>=5)
sub_Sample_2 <- subset(Sample_2,N>=5)
BSobj_sample2_vs_sample1 <- makeBSseqData(list(sub_Sample_1,sub_Sample_2),c("C","N"))
dmltest <- DMLtest(BSobj_sample2_vs_sample1,group1="C",group2="N",smoothing=T)
dmls <- callDML(dmltest, delta=0.2, p.threshold=1e-5)
dmrs <- callDMR(dmltest, p.threshold=1e-5,minlen=50,minCG=3)
write.table(dmls,file=output_DML_file,sep="\t", col.names = T, row.names = F,quote=F)
write.table(dmrs,file=output_DMR_file,sep="\t", col.names = T, row.names = F,quote=F)



tab5rows <- read.table(samp_1_rep1,header=F,sep="\t",nrows=100)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE) #take commandline args
if(length(args) < 5) {
  args <- c("--help")
}



modify_for_DSS <- function(x){
head(tab5rows)
subset_tab5rows <- tab5rows[1,2,3,5,6]
subset_tab5rows <- tab5rows[,c(1,2,3,5,6)]
head(subset_tab5rows)
subset_tab5rows$N <- subset_tab5rows$V5 + subset_tab5rows$V6
head(subset_tab5rows)
subset_again <- subset_tab5rows[,c(1,2,3,4,6)]
head(subset_again)
subset_again <- subset_tab5rows[,c(1,2,3,6,4)]
head(subset_again)
history(max.show=100)
}

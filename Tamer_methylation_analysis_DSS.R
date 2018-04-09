library(DSS)
require(bssseq)

#Reading all the samples and replicate files
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
tab5rows <- read.table(sample_1, header = F,sep="\t", nrows = 100)
classes <- sapply(tab5rows, class)

modify_for_DSS <- function(sample_file,colclasses){
Sample_read <- read.table(sample_file, header = F,sep="\t", colClasses = classes)
subset_sample <- Sample_read[,c(1,2,3,5,6)]
subset_sample$N <- subset_sample$V5 + subset_sample$V6
subset_again <- subset_tab5rows[,c(1,2,6,4)]
colnames(subset_again) <- c("chr","pos","N","X")
sub_Sample_final <- subset(subset_again,N>=5)
return(sub_Sample_final)
}

Sample1_rep1 <- modify_for_DSS(S1_R1,classes)
Sample1_rep2 <- modify_for_DSS(S1_R2,classes)

Sample2_rep1 <- modify_for_DSS(S2_R1,classes)
Sample2_rep2 <- modify_for_DSS(S2_R2,classes)

p.adjusted <- function(pvals,method=c("SLIM","holm","hochberg","hommel",
                                      "bonferroni","BH","BY","fdr","none","qvalue"),
                       n=length(pvals),fdr.level=NULL,pfdr=FALSE,STA=.1,
                       Divi=10,Pz=0.05,B=100,Bplot=FALSE){
  
  method <- match.arg(method)
  
  qvals=switch(method,
               # SLIM function
               SLIM={QValuesfun(pvals,
                                SLIMfunc(pvals,STA=STA,Divi=Divi,Pz=Pz,
                                         B=B,Bplot=Bplot)$pi0_Est)
               },
               # r-base/p-adjust functions
               holm={p.adjust(pvals,method=method,n)},
               hochberg={p.adjust(pvals,method=method,n)},
               hommel={p.adjust(pvals,method=method,n)},
               bonferroni={p.adjust(pvals,method=method,n)},
               BH={p.adjust(pvals,method=method,n)},
               BY={p.adjust(pvals,method=method,n)},
               fdr={p.adjust(pvals,method=method,n)},
               none={p.adjust(pvals,method=method,n)},
               # r-bioconductor/qvalue-package function
               qvalue={qvalue(pvals,fdr.level,pfdr)$qvalues}
  )
}


BSobj_sample2_vs_sample1 <- makeBSseqData(list(Sample1_rep1,Sample1_rep2,Sample2_rep1,Sample2_rep2),c("C1","C2", "N1", "N2"))
dmltest <- DMLtest(BSobj_sample2_vs_sample1,group1=c("C1", "C2"), group2=c("N1", "N2"),smoothing=T)
dmls <- callDML(dmltest, delta=0.2, p.threshold=1e-5)
dmls$qvalue <- p.adjusted(dmls$pvalue,method="bonferroni")
dmrs <- callDMR(dmltest, p.threshold=1e-5,minlen=50,minCG=3)
dmrs$qvalue <- p.adjusted(dmrs$pvalue,method="bonferroni")
write.table(dmls,file=output_DML_file,sep="\t", col.names = T, row.names = F,quote=F)
write.table(dmrs,file=output_DMR_file,sep="\t", col.names = T, row.names = F,quote=F)

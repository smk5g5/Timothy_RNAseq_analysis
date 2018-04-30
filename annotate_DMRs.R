library(annotatr)
all_annots <- builtin_annotations()
allmm10_annots_indices <- grep("mm10",builtin_annotations())
allmm10_annots <- all_annots[allmm10_annots_indices]
sel_mm10_annots <- c("mm10_genes_promoters","mm10_basicgenes","mm10_cpg_islands","mm10_cpg_shelves","mm10_genes_intergenic","mm10_cpg_shores")
mouse_annotations = build_annotations(genome = 'mm10', annotations = sel_mm10_annots)

 # chr      start        end     length        nCG meanMethy1 meanMethy2 diff.Methy 
  #"factor"  "integer"  "integer"  "integer"  "integer"  "numeric"  "numeric"  "numeric" 
  #areaStat 
 #"numeric" 

annotate_DMRs <- function(dm_file,outdmfile){
#chr     start       end length nCG meanMethy1 meanMethy2 diff.Methy  areaStat
extraCols = c(length = 'integer', nCG = 'integer', meanMethy1 = 'numeric',meanMethy2='numeric',diff.Methy="numeric",areaStat="numeric")
pdplus_vs_minus_dml <- read.table(dm_file,header = T,sep = "\t")
myfile<- paste("temp_",dm_file,sep='')
write.table(x = pdplus_vs_minus_dml,file = myfile,sep = "\t",row.names = F,col.names=F,quote = F)
dm_regions = read_regions(con = myfile, genome = 'mm10', extraCols = extraCols, format = 'bed')
dm_annotated = annotate_regions(
    regions = dm_regions,
    annotations = mouse_annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)
write.table(x = df_dm_annotated,file = outdmfile,sep = "\t",row.names = F,quote = F)
if (file.exists(myfile)) file.remove(myfile)
}

            #      chr                    pos                    mu1                    mu2 
             # "factor"              "integer"              "numeric"              "numeric" 
              #    diff                diff.se                   stat                   phi1 
             #"numeric"              "numeric"              "numeric"              "numeric" 
             #     phi2                   pval                    fdr postprob.overThreshold 
             #"numeric"              "numeric"              "numeric"              "numeric"                          
                          
annotate_DMLs <- function(dml_file,outdmlfile){
pdplus_vs_minus_dml <- read.table(dml_file,header = T,sep = "\t")
colnames(pdplus_vs_minus_dml) <- c('chr','pos','mu1','mu2','diff','diff.se','stat','phi1','phi2','pval','fdr','postprob.overThreshold')
pdplus_vs_minus_dml$end <- pdplus_vs_minus_dml$pos + 1
pdplus_vs_minus_dml <- pdplus_vs_minus_dml[c('chr','pos','end','mu1','mu2','diff','diff.se','stat','phi1','phi2','pval','fdr','postprob.overThreshold')]
myfile<- paste("temp_",dml_file,sep='')
write.table(x = pdplus_vs_minus_dml,file = myfile,sep = "\t",row.names = F,col.names=F,quote = F)
extraCols = c(mu1 = 'numeric', mu2 = 'numeric', diff = 'numeric',diff.se='numeric',stat="numeric",phi1="numeric",phi2="numeric",pval="numeric",fdr="numeric",postprob.overThreshold="numeric")
dm_regions = read_regions(con = myfile, genome = 'mm10', extraCols = extraCols, format = 'bed')
dm_annotated = annotate_regions(
    regions = dm_regions,
    annotations = mouse_annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)
write.table(x = df_dm_annotated,file = outdmlfile,sep = "\t",row.names = F,quote = F)
if (file.exists(myfile)) file.remove(myfile)
}
# chr	pos	mu1	mu2	diff	diff.se	stat	phi1	phi2	pval	fdr	postprob.overThreshold

annotate_DMRs("Pdminus_vs_WT_DMR.txt","annotated_Pdminus_vs_WT_DMR.txt")
annotate_DMRs("Pdplus_vs_pdminus_DMR.txt","annotated_Pdplus_vs_pdminus_DMR.txt")
annotate_DMRs("Pdplus_vs_WT_DMR.txt","annotated_Pdplus_vs_WT_DMR.txt")

annotate_DMLs("Pdplus_vs_WT_DML.txt","annotated_Pdplus_vs_WT_DML.txt")
annotate_DMLs("Pdplus_vs_pdminus_DML.txt","annotated_Pdplus_vs_pdminus_DML.txt")
annotate_DMLs("Pdminus_vs_WT_DML.txt","annotated_Pdminus_vs_WT_DML.txt")

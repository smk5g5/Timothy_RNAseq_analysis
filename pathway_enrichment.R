library(gtools)
library(clusterProfiler)
library(org.Mm.eg.db)

pathway_enrichment <- function(dmr_input,pdf_out){
annotated_dmr <- read.table(dmr_input,header = T,sep = "\t")
subset_annot <- annotated_dmr[,c("annot.id","annot.tx_id","annot.gene_id","annot.symbol","meanMethy1","meanMethy2","diff.Methy")]
subse_annot_omit <- na.omit(subset_annot[,c("annot.gene_id","annot.symbol","meanMethy1","meanMethy2","diff.Methy")])
unique_df <- unique(subse_annot_omit)
unique_df$fc <- foldchange(unique_df$meanMethy1,unique_df$meanMethy2)
newdata <- subset(unique_df, fc > 3 | fc < -3)
ggo <- groupGO(gene     = as.character(unique_df$annot.gene_id),
OrgDb    = org.Mm.eg.db,
ont      = "MF",
level    = 12,
readable = TRUE)
pdf(pdf_out)
barplot(ggo, drop=TRUE, showCategory=12)
kk <- enrichKEGG(gene = unique_df,
organism   = 'mmu',
pvalueCutoff = 0.05)
kk <- enrichKEGG(gene = as.character(unique_df$annot.gene_id),
organism   = 'mmu',
pvalueCutoff = 0.05)
barplot(kk)
  dev.off()
return(kk)
}
browseKEGG(kk,'mmu04360')
genelist <- unique_df$fc
names(genelist) <- as.character(unique_df$annot.gene_id)
library("pathview")
source("https://bioconductor.org/biocLite.R")
biocLite("pathview")
library(pathview)
mmu04360 <- pathview(gene.data  = genelist,
pathway.id = "mmu04360",
species    = "mmu",
limit      = list(gene=max(abs(genelist)), cpd=1))
mmu04360
dotplot(kk)








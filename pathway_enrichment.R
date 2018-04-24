library(gtools)
library(clusterProfiler)
library(org.Mm.eg.db)
x <- org.Mm.egCHR
mapped_genes <- mappedkeys(x)
pathway_enrichment <- function(dmr_input,tiffout,rdatafile,sig_path){
annotated_dmr <- read.table(dmr_input,header = T,sep = "\t")
subset_annot <- annotated_dmr[,c("annot.id","annot.tx_id","annot.gene_id","annot.symbol","meanMethy1","meanMethy2","diff.Methy")]
subse_annot_omit <- na.omit(subset_annot[,c("annot.gene_id","meanMethy1","meanMethy2","diff.Methy")])
unique_df <- unique(subse_annot_omit)
aggregate_bygene_df <- aggregate(unique_df,by=list(unique_df$`annot.gene_id`),FUN = "mean")
aggregate_bygene_df$fc <- foldchange(aggregate_bygene_df$meanMethy1,aggregate_bygene_df$meanMethy2)
subset_aggdf <- subset(aggregate_bygene_df[,c("annot.gene_id","meanMethy1","meanMethy2","diff.Methy","fc")], fc > 2 | fc < -2)
ego <- enrichGO(gene          =  as.character(subset_aggdf$annot.gene_id),
                universe      = mapped_genes,
                OrgDb         = org.Mm.eg.db,
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
        readable      = TRUE)  
kk <- enrichKEGG(gene = as.character(subset_aggdf$annot.gene_id),
organism   = 'mmu',
pvalueCutoff = 0.05)
tiff_GO <- unlist(paste("Geneontology_",tiffout,sep=''))
tiff_kegg <- unlist(paste("KEGGenrich_",tiffout,sep=''))
tiff(tiff_GO, width = 20, height = 10, units = 'in', res = 300)
print(barplot(ego, drop=TRUE, showCategory=12))
dev.off()
tiff(tiff_kegg, width = 10, height = 10, units = 'in', res = 300)
print(barplot(kk))
 dev.off()
mykk <- as.data.frame(kk[1:9])
write.table(x = mykk[mykk$Count>0,],file = sig_path,sep = "\t",row.names = F,col.names=T,quote = F)
tiff_map <- unlist(paste("enrichmap_",tiffout,sep=''))
tiff(tiff_map, width = 20, height = 10, units = 'in', res = 300)
print(enrichMap(ego))
dev.off()
save(subset_aggdf,aggregate_bygene_df, ego,kk, file = rdatafile)
return(ego)
}
browseKEGG(kk,'mmu05225')
genelist <- unique_df$fc
names(genelist) <- as.character(unique_df$annot.gene_id)
library(pathview)
mmu04360 <- pathview(gene.data  = genelist,
pathway.id = "mmu04360",
species    = "mmu",
limit      = list(gene=max(abs(genelist)), cpd=1))
mmu04360
dotplot(kk)

#making coldata for deseq2#########
phenotype_data <- read.csv("Mouse_data_timothy.csv",header=F)
colnames(phenotype_data) <- c("sample_name","age_cond")
levels(phenotype_data$age_cond)
#[1] "Old_Control"   "Old_IRP"       "Old_Stretch"   "Young_Control"
#[5] "Young_IRP"     "Young_Stretch"
young_stretch <- phenotype_data[phenotype_data$age_cond == "Young_Stretch",]
old_irp <- phenotype_data[phenotype_data$age_cond == "Old_IRP",]
old_ctrl <- phenotype_data[phenotype_data$age_cond == "Old_Control",]
young_ctrl <- phenotype_data[phenotype_data$age_cond == "Young_Control",]
young_irp <- phenotype_data[phenotype_data$age_cond == "Young_IRP",]
old_irp$age <- "Old"
old_stretch$age <- "Old"
old_ctrl$age <- "Old"
young_ctrl$age <- "Young"
young_irp$age <- "Young"
young_stretch$age <- "Young"
young_ctrl$condition <- "Control"
young_stretch$condition <- "Stretch"
young_irp$condition <- "IRP"
old_ctrl$condition <- "Control"
old_irp$condition <- "IRP"
old_stretch$condition <- "Stretch"
coldata <- rbind(young_ctrl,old_ctrl,young_stretch,old_stretch,young_irp,old_irp)
filelist = list.files(pattern = ".*.count")
count_listdf <- lapply(unlist(filelist),read.table,header=F,sep="\t")
###check if all have equal no of rows###
lapply(count_listdf,nrow)
###extract count columns#######
extracted_counts <- lapply(count_listdf,function(x) x[3])
#######cbind count columns#######
counts_data <- do.call("cbind",extracted_counts)
#######make ensembl gene ids as rownames#######
rownames(counts_data) <- count_listdf[[1]]$V1
colnames(counts_data) <- gsub(".count","",unlist(filelist))
########reorder column names according to coldata#####
reordered_counts <- counts_data[mycoldata$sample_name]




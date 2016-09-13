library(cellrangerRkit)
library(slam)

ts13_1=read.table(file='/data1/ajay/10x/TS13_run1/run1/filtered_2/Matrix_UMI_Gene_Filtered_Clean.tsv',sep='\t',header=T,check.names = F)
ts13_1_filt<-ts13_1
#ts13_1_filt<-ts13_1[rowSums(ts13_1[,2:ncol(ts13_1)])>100,] # Remove genes with less than 100 UMI across all cells
ts13_1_gene_matrix<-ts13_1_filt[,2:ncol(ts13_1_filt)] 
rownames(ts13_1_gene_matrix)<-ts13_1_filt[,1]

#Filter based on gene subset
devpro_genes<-read.table(file='/data1/ajay/10x/TS13_run1/annotations/devpro_ensembl.txt')
ts13_1_gene_matrix_subset<-ts13_1_gene_matrix[,devpro_genes$V1]

ts13_1_gene_matrix_subset<-ts13_1_gene_matrix

ts13_1_gene_matrix_format<-t(ts13_1_gene_matrix_subset)
#ts13_1_gene_matrix_format<-log(ts13_1_gene_matrix_format+0.01,10)
ts13_1_gene_matrix_format_dg<-as(ts13_1_gene_matrix_format,'dgTMatrix')
ts13_1_gene_matrix_format_dg@Dimnames[[1]]<-rownames(ts13_1_gene_matrix_format)
ts13_1_gene_matrix_format_dg@Dimnames[[2]]<-colnames(ts13_1_gene_matrix_format)

gbm <- load_cellranger_matrix(matrix_file_path="/data1/ajay/10x/TS13_run1/matrix.mtx") #load a example gbm object

gbm@mat<-ts13_1_gene_matrix_format_dg # overwrite matrix of gbm object
gbm@gene_symbols<-character(length = 0) 

pca_result <- run_pca(gbm, n_pcs = 10)
#write.table(file = '/data1/ajay/10x/TS13_run1/run1/filtered_2/Matrix_UMI_Gene_Filtered_Clean_tf_10PCA.tsv',pca_result$x,sep='\t',quote = F,col.names = NA)
write.table(file='/data1/ajay/10x/TS13_run1/diff_gene_sets/iterative/10x_run1_ts13_QCd_Wtx_pca_loadings.tsv',pca_result$rotation,sep='\t',quote = F,col.names = NA)
tsne_result <- run_tsne(pca_result,dims=3)
tsne_result_matrix <-tsne_result$Y
rownames(tsne_result_matrix)<-rownames(pca_result$x)
#write.table(file = '/data1/ajay/10x/TS13_run1/run1/filtered_2/Matrix_UMI_Gene_Filtered_Clean_tf_3tSNE.tsv',tsne_result_matrix,sep='\t',quote = F,col.names = NA)

# Perform Mclust on tSNE data
library(mclust)
#mclust analysis
tsne_mclust<-Mclust(tsne_result_matrix,G = 1:150)
summary(tsne_mclust)
plot(tsne_mclust)
mclust_cluster_assignments<-as.data.frame(tsne_mclust$classification)
colnames(mclust_cluster_assignments)<-'mclust=48'
write.csv(mclust_cluster_assignments,file='/data1/ajay/10x/TS13_run1/annotations/TF_tSNE_mclust_assignments.csv',quote = F)

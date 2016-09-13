# Load gene expression matrix #Enter appropriate file name !
ts13_1=read.table(file='/data1/ajay/10x/TS13_run1/matrix.csv',sep=',',header=T,check.names = F)
# Load QCed barcodes
cell_qc=read.table(file='/data1/ajay/10x/TS13_run1/run1/filtered_2/Matrix_UMI_Gene_Filtered_Clean.tsv',sep='\t',header=T,check.names = F)

#Filter based on gene subset
devpro_genes<-read.table(file='/data1/ajay/10x/TS13_run1/annotations/devpro_ensembl.txt',na.strings=c("", "NA"))
devpro_genes<-as.data.frame(unique(devpro_genes$V1))


ts13_1_filt<-merge(ts13_1,devpro_genes,by.x=1,by.y=1)
#ts13_1_filt<-ts13_1[rowSums(ts13_1[,3:ncol(ts13_1)])>0,] # Remove genes with no expression across all samples
ts13_1_gene_matrix<-ts13_1_filt[,3:ncol(ts13_1_filt)] 
rownames(ts13_1_gene_matrix)<-ts13_1_filt[,1]
ts13_1_gene_matrix2<-ts13_1_gene_matrix[,as.character(cell_qc$BarCode)]

ensembl_gene_name_df<-as.data.frame(ts13_1_filt[,2])
rownames(ensembl_gene_name_df)<-ts13_1_filt[,1]
colnames(ensembl_gene_name_df)<-'gene_name'

#Load Annotation table

ts13_1_annot=read.table(file='/data1/ajay/10x/TS13_run1/annotations/Matrix_UMI_Gene_Filtered_Clean_devpro-genes_10PCA_Clustered.csv',sep=',',header=T,check.names = F,row.names = 1)
#colnames(ts13_1_annot)<-'W=20'

ts13_1_gene_matrix_subset<-ts13_1_gene_matrix2[,rownames(ts13_1_annot)]

if (!identical(rownames(ts13_1_annot),colnames(ts13_1_gene_matrix_subset))){  #check if rownames in annotation matrix are in the same order
  ts13_1_annot<-ts13_1_annot[match(colnames(ts13_1_gene_matrix_subset),as.character(ts13_1_annot$`Bar Code`)),] #if not in same order, reorder matrix based on order from gene expresion matrix
}


for (numk in c('K20','K30','K36',19:24))   #hard code the columns that need to iterated
 {
#numk='mclust=48'
  df_umi_count_avg<-list()
  colnames_vec<-c()
  n=1
for (x in sort(unique(ts13_1_annot[,numk]))){
  #for (shuffle in 1:10){   comment out for bootstrap
    sample_names<-rownames(ts13_1_annot)[ts13_1_annot[,numk]==x]
    #sample_names<-sample(sample_names,ceiling(length(sample_names)*0.75)) # comment out for bootstrap
    df_umi_count_avg[[n]]<-apply(ts13_1_gene_matrix_subset[,sample_names],1,mean)
    
    #colnames_vec<-append(colnames_vec,paste('cluster_',x,'_boot_',shuffle,sep='')) # comment out for bootstrap
    colnames_vec<-append(colnames_vec,paste('cluster_',x,sep='')) # comment in for bootstrap
    n=n+1
  #}
}
  df_umi_count_avg<-as.data.frame(df_umi_count_avg)
  colnames(df_umi_count_avg)<-colnames_vec
  df_umi_count_avg$gene_name<-ensembl_gene_name_df$gene_name
  write.table(file=paste('/data1/ajay/10x/TS13_run1/diff_gene_sets/pca/10x_run1_ts13_QCd_devpro_pca_',numk,'_cluster_avg.csv',sep=''),df_umi_count_avg,sep=',',quote = F,row.names = T,col.names = NA)
}

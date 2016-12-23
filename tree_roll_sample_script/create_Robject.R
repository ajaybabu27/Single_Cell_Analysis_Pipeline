#!/usr/bin/env Rscript

#Author: Ajay
#Description: Generate background distributions from single cell gene expression data. Measure euclidean distances (to Anatomical Structures) for all
#background vectors. The output R object (data frame) is used as background model to infer p-values when used in Tree roll matching. 

library(vegan)
library('hopach')

args <- commandArgs(TRUE)




cell_file=args[1]
file=args[2]

term_data=read.table(file, header = TRUE, sep="\t",row.names=1)

bootstraps=10

cell_data_back=read.csv(file = cell_file,row.names = 1)
cell_data_back=cell_data_back[cell_data_back$QC=='live',]
cell_data_back<-cell_data_back[,1:96]
cell_random<- permatfull(cell_data_back, fixedmar="columns", times = bootstraps, mtype = "prab")

cell_match_data_frame_back_all=data.frame(matrix(NA,ncol=nrow(term_data),nrow=0))

colnames(cell_match_data_frame_back_all)<-row.names(term_data)
for (randomize in 1:bootstraps){
  print(randomize)
  cell_match_data_frame_back=matrix(0,ncol=nrow(term_data),nrow=nrow(cell_data_back))
  colnames(cell_match_data_frame_back)=row.names(term_data)
  rownames(cell_match_data_frame_back)=rownames(cell_data_back)
  for (dist_meth in c('euclid'))
    
  {
    for(n in 1:nrow(cell_data_back))
    {
      vector=as.vector(cell_random$perm[[randomize]][n,],mode='numeric')
      #cell_dist<-distancevector(term_data[1:ncol(term_data)-1],vector, d=dist_meth, na.rm=TRUE)
      cell_dist<-distancevector(term_data[,1:96],vector, d=dist_meth, na.rm=TRUE)
      cell_match_data_frame_back[n,]<-cell_dist
    }
  }    
  
  cell_match_data_frame_back<-as.data.frame(cell_match_data_frame_back)
  cell_match_data_frame_back_all<-rbind(cell_match_data_frame_back_all,cell_match_data_frame_back)  
}

back_data<-cell_match_data_frame_back_all
save(back_data,file='/home/single_cell_dev/lineage_map/tree_roll_sample_script/output/back.rds')


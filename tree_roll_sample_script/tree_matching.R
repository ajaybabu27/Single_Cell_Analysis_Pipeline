#!/usr/bin/env Rscript

#Author: Ajay
#Description: Computes the p-value for each single cell & anatomical strucuture pair at a given tree level. 
#Finally the anatomical structure with the lowest p-value is chosen.


args <- commandArgs(TRUE)


load(args[1]) # Load Background model R Object

f <- function(x,y){
  #test <- t.test(x,y,alternative = 'less')
  pval_vec=c()
  a <- mean(y)
  s <- sd(y)
  n <- length(y)
  
  for (n in 1:length(x)){
    
    xbar <- x[n]
    pval=2*(1-pnorm(xbar,mean=a,sd=s/sqrt(n)))  
    #pval=sum(y <= x[n])/(nrow(back_data))
    pval_vec<-append(pval_vec,pval)
  }
  data.frame(pval_vec)
  #data.frame(pval = test$p.value)
  
}

library('hopach')

ref_file=args[2] # Load temp file containing single cell and anatomical structure nodes being compared against.
term_data=read.table(ref_file, header = TRUE, sep="\t",row.names=1)
cell_data=term_data[1,]
ref_data=term_data[2:nrow(term_data),]

cell_match_data_frame=matrix(0,ncol=nrow(cell_data),nrow=nrow(ref_data))

for(n in 1:nrow(cell_data))
{
  vector=as.vector(cell_data[n,],mode='numeric')
  cell_dist<-distancevector(ref_data,vector, d='euclid', na.rm=TRUE)
  cell_match_data_frame[,n]<-cell_dist
}

#cell_match_data_frame<-rescale(t(cell_match_data_frame),c(0.1,1))
cell_match_data_frame=t(cell_match_data_frame)
colnames(cell_match_data_frame)<-rownames(ref_data)

p_vals<-as.matrix(as.data.frame(sapply(intersect(colnames(cell_match_data_frame),colnames(back_data)), function(x) f(cell_match_data_frame[,x], back_data[,x]))))



p_vals_fdr<-apply(p_vals,1, function(x) p.adjust(x, method='fdr',n=length(x))) #Need to work on this !
inds = which(p_vals_fdr == min(p_vals_fdr), arr.ind=TRUE)
rnames=rownames(p_vals_fdr)[inds[,1]]
print(strsplit(rnames,'.p')[[1]][1])



library("R.matlab")
library("readr")

#######define function###################
#calculate gene_wise correlation
get.cor.gene <- function(X, Y) {
  sapply(1:nrow(X), function(i) cor(X[i, ], Y[i, ]))
}

#calculate cell_wise correlation
get.cor.cell <- function(X, Y) {
  sapply(1:ncol(X), function(i) cor(X[, i], Y[, i]))
}

## Correlation matrix distance(CMD)
calc_cmd <- function(R1, R2) {
  traceR1R2 <- sum(diag(crossprod(R1, R2)))
  R1.norm <- norm(R1, type = "F")
  R2.norm <- norm(R2, type = "F")
  return(1-traceR1R2/(R1.norm*R2.norm))
}

###############################################################################
t1=proc.time()
method = c('dropout','Bubble','SAVER','ALRA','VIPER','scimpute','bayNorm','auto','DCA','scrabble')
dataset = c('sim_data1.1','sim_data1.2','sim_data2.1','sim_data2.2','sim_data3.1','sim_data3.2','sim_data4.1','sim_data4.2','sim_data5.1','sim_data5.2')
#dataset = c('sim_data1','sim_data2','sim_data3','sim_data4','sim_data5')
method_number=length(method)
dataset_number=length(dataset)
gene_cor_list=list()
cell_cor_list=list()
gene_cmd_list=list()
cell_cmd_list=list()
for (ii in 1:dataset_number){
  print(dataset[ii])
  data_ref_dir=paste('./dataset/',dataset[ii],'_true.txt',sep='')
  data_ref=read.csv(data_ref_dir, header = T, row.names = 1,sep='\t')
  data_ref=t(as.matrix(data_ref))
  for (i in 1:method_number){
    print(method[i])
    if (method[i] == 'dropout'){
      data_imp_dir = paste('./dataset/',dataset[ii],'_row_col.txt',sep='')
      data_imp=read.csv(data_imp_dir, header = T, row.names = 1,sep='\t')
      data_imp=as.matrix(data_imp)
    }else if(method[i] == 'Bubble'){
      data_imp_dir = paste('./imputed_data/',method[i],'/',dataset[ii],'.txt',sep='')
      data_imp=read.csv(data_imp_dir, header = T, row.names = 1,sep='\t')
      data_imp=as.matrix(data_imp)
    }else if(method[i] == 'scimpute'){
      data_imp_dir = paste('./imputed_data/',method[i],'/',dataset[ii],'.txt',sep='')
      data_imp=read.csv(data_imp_dir,header = T, row.names = 1,sep=' ')
      data_imp=as.matrix(data_imp)
    }else if(method[i] == 'DCA'){
      data_imp_dir = paste('./imputed_data/',method[i],'/',dataset[ii],'_results','/','mean.tsv',sep='')
      data_imp=read_tsv(data_imp_dir, col_names = TRUE)
      data_imp=as.data.frame(data_imp)
      data_imp=as.matrix(data_imp)
      data_imp=data_imp[-1,]
      gene_col=data_imp[,1]
      gene_rowname=colnames(data_ref)[1:1500]
      gene_rowname=as.numeric(gene_rowname)
      miss_geneindex=setdiff(gene_rowname,gene_col)
      data_ref=data_ref[-miss_geneindex,]
      data_imp=data_imp[,-1]
    }else if(method[i] == 'auto'){
      data_imp_dir = paste('./imputed_data/',method[i],'/',dataset[ii],'.mat',sep='')
      data_imp=readMat(data_imp_dir)
      col_num=dim(data_imp$arr)[3]
      data_imp=matrix(aperm(data_imp$arr, c(1, 3, 2)), ncol = col_num)
      data_imp=t(data_imp)
    }else if(method[i] == 'scrabble'){
      data_imp_dir = paste('./imputed_data/',method[i],'/',dataset[ii],'.txt',sep='')
      data_imp=read.csv(data_imp_dir, header = T, row.names = 1,sep=' ')
      data_imp=as.matrix(data_imp)
    }else{
      data_imp_dir = paste('./imputed_data/',method[i],'/',dataset[ii],'.txt',sep='')
      data_imp=read.csv(data_imp_dir, header = T, row.names = 1,sep='\t')
      data_imp=as.matrix(data_imp)
    }
    
    print(paste("imp_data_dim:",dim(data_imp)))
    print(paste("ref_data_dim:",dim(data_ref)))
    
    #calculate gene_wise correlation
    gene_cor<-get.cor.gene(data_ref, data_imp)
    gene_cor_mean<-mean(gene_cor,na.rm=TRUE)
    gene_cor_list[[method[i]]]<-append(gene_cor_list[[method[i]]],gene_cor_mean)
    #paste("gene_cor:",gene_cor)
    print(paste("gene_cor_mean:",gene_cor_mean))

    #calculate cell_wise correlation
    cell_cor<-get.cor.cell(data_ref, data_imp)
    cell_cor_mean<-mean(cell_cor,na.rm=TRUE)
    cell_cor_list[[method[i]]]<-append(cell_cor_list[[method[i]]],cell_cor_mean)
    #paste("cell_cor:",cell_cor)
    print(paste("cell_cor_mean:",cell_cor_mean))

    #calculate gene correlation matrix(cor())

    zero_rowindex1=which(apply(data_imp,1,sd)==0)
    zero_rowindex2=which(apply(data_ref,1,sd)==0)
    gene_index=union(zero_rowindex1,zero_rowindex2)
    if(length(gene_index)>0){
      data_imp1=data_imp[-gene_index,]
      data_ref1=data_ref[-gene_index,]
      print(dim(data_imp1))
      print(dim(data_ref1))
    }else{
      data_imp1=data_imp
      data_ref1=data_ref
      print(dim(data_imp1))
      print(dim(data_ref1))
    }
    R1=cor(t(data_ref1))
    R2=cor(t(data_imp1))
    gene_cmd=calc_cmd(R1,R2)
    gene_cmd_list[[method[i]]]<-append(gene_cmd_list[[method[i]]],gene_cmd)
    #paste('gene_cmd:',gene_cmd)
    print(paste("gene_cmd:",gene_cmd))

    #calculate cell correlation matrix(cor())
    zero_colindex1=which(apply(data_imp,2,sd)==0)
    zero_colindex2=which(apply(data_ref,2,sd)==0)
    cell_index=union(zero_colindex1,zero_colindex2)
    if(length(cell_index)>0){
      data_imp1=data_imp[,-col_index]
      data_ref1=data_ref[,-col_index]
      print(dim(data_imp1))
      print(dim(data_ref1))
    }else{
      data_imp1=data_imp
      data_ref1=data_ref
      print(dim(data_imp1))
      print(dim(data_ref1))
    }
    R1=cor(data_ref1)
    R2=cor(data_imp1)
    cell_cmd=calc_cmd(R1,R2)
    cell_cmd_list[[method[i]]]<-append(cell_cmd_list[[method[i]]],cell_cmd)
    #paste('cell_cmd:',cell_cmd)
    print(paste("cell_cmd:",cell_cmd))
  }
}

###save datasets####
write.table(gene_cor_list,paste0("./Cmd-Cor/","gene_cor.txt",sep=''))
write.table(cell_cor_list,paste0("./Cmd-Cor/","cell_cor.txt",sep=''))
write.table(gene_cmd_list,paste0("./Cmd-Cor/","gene_cmd.txt",sep=''))
write.table(cell_cmd_list,paste0("./Cmd-Cor/","cell_cmd.txt",sep=''))
print('Finish saving')
t2=proc.time()
t=t2-t1
print(paste0('running time：',t[3][[1]],'second'))


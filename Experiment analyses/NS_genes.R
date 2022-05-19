library(readr)
library(R.matlab) 

###############################################################################
t1=proc.time()
method = c('dropout','Bubble','SAVER','scimpute','VIPER','bayNorm','ALRA','auto',"scrabble",'DCA' )
dataset = c('sim_data1','sim_data2','sim_data3','sim_data4','sim_data5','sim_data1.1','sim_data1.2','sim_data2.1','sim_data2.2','sim_data3.1','sim_data3.2','sim_data4.1','sim_data4.2','sim_data5.1','sim_data5.2')
method_number=length(method)
dataset_number=length(dataset)
for (ii in 1:dataset_number){
  NS_GENE_LIST=list()
  print(dataset[ii])

  #####sim_data1-sim_data5.2##
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
      rownames(data_ref)=c(0:1499)
      data_imp=data_imp[-1,]
      gene_col=data_imp[,1]
      gene_rowname=rownames(data_ref)[1:1500]
      gene_rowname=as.numeric(gene_rowname)
      intersect_geneindex=intersect(gene_rowname,gene_col)
      intersect_geneindex=as.character(intersect_geneindex)
      print(length(intersect_geneindex))
      data_ref=data_ref[intersect_geneindex,]
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
    # ##############wilcoxon test################
    # ## wilcoxon
    pval = sapply(1:nrow(data_imp), function(i1) {
      gene_one = data_ref[i1,]
      gene_two = data_imp[i1,]
      if (identical(gene_one,gene_two)) {
        1
      }else{
        res = wilcox.test(gene_one,gene_two)
        res$p.value
      }
    })
    # ##output the number of diffgene that p.adjust less than 0.01
    fdr = p.adjust(pval, method='fdr')
    fdr_list=list(fdr)
    diffgene_raw=which(fdr_list[[1]] > 0.01)
    NS_GENE_NUM=length(diffgene_raw)
    print(paste0('NS_gene_num：',NS_GENE_NUM))
    NS_GENE_LIST=append(NS_GENE_LIST,NS_GENE_NUM)
  }
  print(NS_GENE_LIST)
  write.table(NS_GENE_LIST,paste0("./NS_GENE/",dataset[ii],"_NSgene_num.txt",sep=''))
}
t2=proc.time()
t=t2-t1
print(paste0('NS_gene_num_test time：',t[3][[1]],'second'))

library(readr)
t1=proc.time()
method = c('dropout','Bubble','SAVER','ALRA','VIPER','scimpute','bayNorm','auto','DCA',"scrabble")
dataset = c('sim_data5.2')
dataset_number=length(dataset)
method_number=length(method)
####Initial list####
diffgene_index_list=list(list())
diffgene_num_list=list(list())
diffgene_overlap_list=list(list())
f1_score_list=list(list())

if (method_number>length(diffgene_index_list) | method_number>length(diffgene_num_list) | method_number>length(diffgene_overlap_list) | method_number>length(f1_score_list))
  stop("there is a mistake,the number of methods not equal the length of lists")
for (ii in 1:dataset_number){
  print(dataset[ii])
  # Load the single-cell label
  label_dir = paste('./label/',dataset[ii],'_label.txt',sep='')
  labelName = read.table(label_dir,sep = "\t",header = FALSE)
  print(dim(labelName))
  labelName=as.numeric(labelName)
  celltype = unique(labelName)
  celltype_num = length(celltype)  
  print(celltype_num)
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
      gene_col=data_imp[,1]
      rownames(data_imp)=gene_col
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
    print(paste("data_imp_dim:",dim(data_imp)))
    print(paste("data_ref_dim:",dim(data_ref)))
    ####test all cellline in the dataset
    for (j in 1:(celltype_num-1)) {
      for (jj in (j+1):celltype_num){
        type=paste('celltype',celltype[j],' vs ','celltype',celltype[jj],sep='')
        print(type)
        ## wilcoxon
        pval = sapply(1:nrow(data_imp), function(i1) {
          celltype_one = data_imp[i1,labelName %in% celltype[j]]
          celltype_two = data_imp[i1,labelName %in% celltype[jj]]
          names(celltype_one) <- names(celltype_two) <- NULL
          if (identical(celltype_one,celltype_two)) {
            1
          }else{
            res = wilcox.test(celltype_one,celltype_two)
            res$p.value
          }
        })
        ##output the number of diffgene that p.adjust less than 0.05
        fdr = p.adjust(pval, method='fdr')
        fdr_list=list(fdr)
        diffgene=which(fdr_list[[1]] < 0.05)
        diffgene_index_list[[i]][[type]]<-append(diffgene_index_list[[i]][[type]],diffgene)
        diffgene_num=length(diffgene)
        print(diffgene_num)
        diffgene_num_list[[i]][[type]]<-append(diffgene_num_list[[i]][[type]],diffgene_num)
        diffgene_gold=diffgene_gold_list[[type]]
        diffgene_overlap=length(intersect(diffgene_gold,diffgene))
        print(diffgene_overlap)
        diffgene_overlap_list[[i]][[type]]<-append(diffgene_overlap_list[[i]][[type]],diffgene_overlap)
        precison=diffgene_overlap/diffgene_num
        recall=diffgene_overlap/length(diffgene_gold)
        f1_score=(2*precison*recall)/(precison+recall)
        print(f1_score)
        f1_score_list[[i]][[type]]<-append(f1_score_list[[i]][[type]],f1_score)
      }
    }
  }
}
t2=proc.time()
t=t2-t1
print(paste0('diffgene_test：',t[3][[1]],'second'))
###save datasets####
write.table(f1_score_list,paste0("./DEG-test/",method[1],'_',dataset[1],"_F1.txt",sep=''))
write.table(diffgene_num_list,paste0("./DEG-test/",method[1],'_',dataset[1],"_diffgene_num.txt",sep=''))
write.table(diffgene_overlap_list,paste0("./DEG-test/",method[1],'_',dataset[1],"_diffgene_overlap.txt",sep=''))
print(paste0(dataset[1],'  Finish saving'))











library(stats)
t1=proc.time()
dataset = c('sim_data5.2')
dataset_number=length(dataset)
#Initial list
diffgene_gold_list=list()
diffgene_gold_num_list=list()
for (ii in 1:dataset_number){
  print(dataset[ii])
  # Load the true_data label
  label_dir = paste('./label/',dataset[ii],'_label.txt',sep='')
  labelName = read.table(label_dir,sep = "\t",header = FALSE)
  print(dim(labelName))
  labelName=as.numeric(labelName)
  celltype = unique(labelName)
  celltype_num = length(celltype)  
  print(celltype_num)
  ####Load the true_data####
  data_ref_dir=paste('./dataset/',dataset[ii],'_true.txt',sep='')
  data_ref=read.csv(data_ref_dir, header = T, row.names = 1,sep='\t')
  data_ref=t(as.matrix(data_ref))
  print(dim(data_ref))
  ####test all cellline in the dataset
  for (j in 1:(celltype_num-1)) {
    for (jj in (j+1):celltype_num){
      print(paste(celltype[j],'vs',celltype[jj]))
      ## wilcoxon
      pval = sapply(1:nrow(data_ref), function(i1) {
        celltype_one = data_ref[i1,labelName %in% celltype[j]]
        celltype_two = data_ref[i1,labelName %in% celltype[jj]]
        names(celltype_one) <- names(celltype_two) <- NULL
        if (identical(celltype_one,celltype_two)) {
          1
        }else{
          res = wilcox.test(celltype_one,celltype_two)
          res$p.value
        }
      })
      ##output the number of diffgene that p.adjust less than 0.05
      fdr_gold = p.adjust(pval, method='fdr')
      fdr_list_gold=list(fdr_gold)
      diffgene_gold_index=which(fdr_list_gold[[1]] < 0.05)
      diffgene_gold_list[[paste('celltype',celltype[j],' vs ','celltype',celltype[jj],sep='')]]=append(diffgene_gold_list[[paste('celltype',celltype[j],' vs ','celltype',celltype[jj],sep='')]],diffgene_gold_index)
      diffgene_gold_num=length(diffgene_gold_index)
      print(diffgene_gold_num)
      diffgene_gold_num_list[[paste('celltype',celltype[j],' vs ','celltype',celltype[jj],sep='')]]=append(diffgene_gold_num_list[[paste('celltype',celltype[j],' vs ','celltype',celltype[jj],sep='')]],diffgene_gold_num)
    }
  }
}
###save datasets####
write.table(diffgene_gold_num_list,paste0("./DEG-test/",dataset[ii],"_diffgene_gold_num.txt",sep=''))
print(paste0(dataset[ii],'  Finish saving'))
t2=proc.time()
t=t2-t1
print(paste0('diffgene_gold_test：',t[3][[1]],'second'))


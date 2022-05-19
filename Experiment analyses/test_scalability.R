library(ggplot2)
method = c('autoimpute','SAVER','scImpute','VIPER','bayNorm','ALRA','DCA','AutoImpute','SCRABBLE')
method_number=length(method)
cell_num= c(1000,5000,50000)
cell_num<- log10(cell_num)
time_data=read.table("./time-memory/time.csv",sep=",", header = TRUE, row.names = 1)
for (i in 1:method_number){
  myfit<-lm(as.numeric(time_data[i,])~cell_num)
  scalability=summary(myfit)$coefficients[2,1]
  print(paste0(method[i],' ',scalability))
}

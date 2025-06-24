source("2d.R")
for(i in 1:50){
  set.seed(i+10086)
  #training
  n0=n1=n2=100; M1=M2=10; eps=0.5; R=5; n.prop=0.5*min(c(n0,n1,n2))
  Data=Data_1(n0, n1, n2, M1, M2, eps, R)
  J1=J2=10; J.opt=10
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  Score=lapply(Data[1:3], score_extract, J1, J2, J.opt, M1, M2)
  label=Data$noise_label
  Score=do.call(rbind, Score)
  for(j in 1:R){
    ind=c(sample(1:n0,n.prop,replace = FALSE), sample((n0+1):(n0+n1),n.prop,replace = FALSE), sample((n0+n1+1):(n0+n1+n2),n.prop,replace = FALSE))
    label[-ind,j] = NaN
  }
  name.1=paste0("trainX",i,".csv")
  name.2=paste0("trainY",i,".csv")
  write.csv(Score, file=name.1, row.names = FALSE, col.names = FALSE)
  write.csv(label, file=name.2, row.names = FALSE, col.names = FALSE)
  #testing
  n0=n1=n2=30; n.prop=0.5*min(c(n0,n1,n2))
  Data=Data_1(n0, n1, n2, M1, M2, eps, R)
  Score=lapply(Data[1:3], score_extract, J1, J2, J.opt, M1, M2)
  label=Data$noise_label
  Score=do.call(rbind, Score)
  for(j in 1:R){
    ind=c(sample(1:n0,n.prop,replace = FALSE), sample((n0+1):(n0+n1),n.prop,replace = FALSE), sample((n0+n1+1):(n0+n1+n2),n.prop,replace = FALSE))
    label[-ind,j] = NaN
  }
  name.1=paste0("testX",i,".csv")
  name.2=paste0("testY",i,".csv")
  write.csv(Score, file=name.1, row.names = FALSE, col.names = FALSE)
  write.csv(label, file=name.2, row.names = FALSE, col.names = FALSE)
}

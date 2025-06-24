#######Scenario 1 for 2D##########
Data_1=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(8, 7, 6, 5, 4); eigen1=c(5, 4, 3, 2, 1); eigen2=c(2.5, 2, 1.5, 1, .5)
  mu0=c(4,4,3,3,3); mu1=c(-1,-1,-1,-1,-1); mu2=c(0,0,0,0,0)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rnorm(n1, mu1[1], eigen1[1]), rnorm(n1, mu1[2], eigen1[2]), rnorm(n1, mu1[3], eigen1[3]), rnorm(n1, mu1[4], eigen1[4]), rnorm(n1, mu1[5], eigen1[5]))
  xi2=cbind(rnorm(n2, mu2[1], eigen2[1]), rnorm(n2, mu2[2], eigen2[2]), rnorm(n2, mu2[3], eigen2[3]), rnorm(n2, mu2[4], eigen2[4]), rnorm(n2, mu2[5], eigen2[5]))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(NA,3,3)
  
  diag(CM)=1-eps
  CM[row(CM) != col(CM)] <- eps/(3-1)
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}


Fourier=function(s, M, j){
  k=j %/% 2
  
  if(j==1){
    return(1)
  }else if(j %% 2 == 0){
    return(sqrt(2/M)*cos(2*pi*k*s))
  }else if(j %% 2 != 0){
    return(sqrt(2/M)*sin(2*pi*k*s))
  }
}



score_extract=function(Data, J1, J2, J.opt, M1, M2){

  J=J1*J2
  
  phi1=matrix(NA, M1, J1)
  for(m in 1:M1){
    for(j in 1:J1){
      phi1[m, j]=Fourier(S1[m], M1, j)
    }
  }
  phi2=matrix(NA, M2, J2)
  for(m in 1:M2){
    for(j in 1:J2){
      phi2[m, j]=Fourier(S2[m], M2, j)
    }
  }
  
  phi=kronecker(t(phi2), t(phi1))
  
  
  Score=matrix(NA, nrow(Data), J)
  
  for(i in 1:nrow(Data)){
    for(j in 1:J){
      Score[i, j] = mean(Data[i,]*phi[j,])
    }
  }
  
  
  Score=Score[,1:J.opt]
  
  return(Score=Score)
  
} 
  

#######Scenario 2 for 2D##########
Data_2=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(5,4,3,2,1); eigen1=c(5, 4, 3, 2, 1)/2; df2=c(3,5,7,9,11)
  mu0=c(-1, -1, -1, -1, -1); mu1=c(0,0,0,0,0)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rnorm(n1, mu1[1], eigen1[1]), rnorm(n1, mu1[2], eigen1[2]), rnorm(n1, mu1[3], eigen1[3]), rnorm(n1, mu1[4], eigen1[4]), rnorm(n1, mu1[5], eigen1[5]))
  xi2=cbind(rt(n2, df2[1],0), rnorm(n2, df2[2],0), rnorm(n2, df2[3],0), rnorm(n2, df2[4],0), rnorm(n2, df2[5],0))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(NA,3,3)
  
  diag(CM)=1-eps
  CM[row(CM) != col(CM)] <- eps/(3-1)
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}








#######Scenario 3 for 2D##########
Data_3=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(5,4,3,2,1)/2; df1=c(2,3,4,5,6); df2=c(3,5,7,9,11)
  mu0=c(0,0,0,0,0); ncp1=c(1,1,1,1,1); ncp2=c(3,3,3,3,3)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rt(n1, df1[1],ncp1[1]), rnorm(n1, df1[2],ncp1[2]), rnorm(n1, df1[3],ncp1[3]), rnorm(n1, df1[4],ncp1[4]), rnorm(n1, df1[5], ncp1[5]))
  xi2=cbind(rt(n2, df2[1],ncp2[1]), rnorm(n2, df2[2],ncp2[2]), rnorm(n2, df2[3],ncp2[3]), rnorm(n2, df2[4],ncp2[4]), rnorm(n2, df2[5], ncp2[5]))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(NA,3,3)
  
  diag(CM)=1-eps
  CM[row(CM) != col(CM)] <- eps/(3-1)
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}









#######Scenario 4 for 2D##########
Data_4=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(5,4,3,2,1)/2; r1=c(0.1,0.3,0.5,0.7,0.9); df2=c(3,5,7,9,11)
  mu0=c(0,0,0,0,0); ncp2=c(3,3,3,3,3)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rexp(n1, r1[1]), rexp(n1, r1[2]), rexp(n1, r1[3]), rexp(n1, r1[4]), rexp(n1, r1[5]))
  xi2=cbind(rt(n2, df2[1],ncp2[1]), rnorm(n2, df2[2],ncp2[2]), rnorm(n2, df2[3],ncp2[3]), rnorm(n2, df2[4],ncp2[4]), rnorm(n2, df2[5], ncp2[5]))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(NA,3,3)
  
  diag(CM)=1-eps
  CM[row(CM) != col(CM)] <- eps/(3-1)
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}






#######Scenario 5 for 2D##########
Data_5=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(8, 7, 6, 5, 4); eigen1=c(5, 4, 3, 2, 1); eigen2=c(2.5, 2, 1.5, 1, .5)
  mu0=c(4,4,3,3,3); mu1=c(-1,-1,-1,-1,-1); mu2=c(0,0,0,0,0)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rnorm(n1, mu1[1], eigen1[1]), rnorm(n1, mu1[2], eigen1[2]), rnorm(n1, mu1[3], eigen1[3]), rnorm(n1, mu1[4], eigen1[4]), rnorm(n1, mu1[5], eigen1[5]))
  xi2=cbind(rnorm(n2, mu2[1], eigen2[1]), rnorm(n2, mu2[2], eigen2[2]), rnorm(n2, mu2[3], eigen2[3]), rnorm(n2, mu2[4], eigen2[4]), rnorm(n2, mu2[5], eigen2[5]))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(0,3,3)
  
  diag(CM)=1-eps
  CM[1,2] <- eps
  CM[2,3] <- eps
  CM[3,1] <- eps
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}


#######Scenario 6 for 2D##########
Data_6=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(5,4,3,2,1); eigen1=c(5, 4, 3, 2, 1)/2; df2=c(3,5,7,9,11)
  mu0=c(-1, -1, -1, -1, -1); mu1=c(0,0,0,0,0)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rnorm(n1, mu1[1], eigen1[1]), rnorm(n1, mu1[2], eigen1[2]), rnorm(n1, mu1[3], eigen1[3]), rnorm(n1, mu1[4], eigen1[4]), rnorm(n1, mu1[5], eigen1[5]))
  xi2=cbind(rt(n2, df2[1],0), rnorm(n2, df2[2],0), rnorm(n2, df2[3],0), rnorm(n2, df2[4],0), rnorm(n2, df2[5],0))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(0,3,3)
  
  diag(CM)=1-eps
  CM[1,2] <- eps
  CM[2,3] <- eps
  CM[3,1] <- eps
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}








#######Scenario 7 for 2D##########
Data_7=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(5,4,3,2,1)/2; df1=c(2,3,4,5,6); df2=c(3,5,7,9,11)
  mu0=c(0,0,0,0,0); ncp1=c(1,1,1,1,1); ncp2=c(3,3,3,3,3)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rt(n1, df1[1],ncp1[1]), rnorm(n1, df1[2],ncp1[2]), rnorm(n1, df1[3],ncp1[3]), rnorm(n1, df1[4],ncp1[4]), rnorm(n1, df1[5], ncp1[5]))
  xi2=cbind(rt(n2, df2[1],ncp2[1]), rnorm(n2, df2[2],ncp2[2]), rnorm(n2, df2[3],ncp2[3]), rnorm(n2, df2[4],ncp2[4]), rnorm(n2, df2[5], ncp2[5]))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(0,3,3)
  
  diag(CM)=1-eps
  CM[1,2] <- eps
  CM[2,3] <- eps
  CM[3,1] <- eps
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}









#######Scenario 8 for 2D##########
Data_8=function(n0, n1, n2, M1, M2, eps, R){
  S1=seq(0, 1, length.out=M1)
  S2=seq(0, 1, length.out=M2)
  
  
  
  eigen0=c(5,4,3,2,1)/2; r1=c(0.1,0.3,0.5,0.7,0.9); df2=c(3,5,7,9,11)
  mu0=c(0,0,0,0,0); ncp2=c(3,3,3,3,3)
  
  xi0=cbind(rnorm(n0, mu0[1], eigen0[1]), rnorm(n0, mu0[2], eigen0[2]), rnorm(n0, mu0[3], eigen0[3]), rnorm(n0, mu0[4], eigen0[4]), rnorm(n0, mu0[5], eigen0[5]))
  xi1=cbind(rexp(n1, r1[1]), rexp(n1, r1[2]), rexp(n1, r1[3]), rexp(n1, r1[4]), rexp(n1, r1[5]))
  xi2=cbind(rt(n2, df2[1],ncp2[1]), rnorm(n2, df2[2],ncp2[2]), rnorm(n2, df2[3],ncp2[3]), rnorm(n2, df2[4],ncp2[4]), rnorm(n2, df2[5], ncp2[5]))
  
  
  SS1=rep(S1, M2)
  SS2=rep(S2, each=M1)
  
  
  BB1=SS1; BB2=SS2; 
  
  BB3=(SS1)*(SS2); BB4=(SS1)^2; BB5=(SS2)^2 
  
  #BB6=(SS1)^2*(SS2); BB7=(SS1)*(SS2)^2; BB8=(SS1)^3; BB9=(SS2)^3
  
  
  BB=rbind(BB1, BB2, BB3, BB4, BB5)
  
  X_0=as.matrix(xi0%*%BB)
  X_1=as.matrix(xi1%*%BB)
  X_2=as.matrix(xi2%*%BB)
  
  
  CM=matrix(0,3,3)
  
  diag(CM)=1-eps
  CM[1,2] <- eps
  CM[2,3] <- eps
  CM[3,1] <- eps
  
  
  noise_label=matrix(NA, n0+n1+n2, R)
  for(r in 1:R){
    for(k in 1:3){
      noise_label[,r]=c(sample(c(0,1,2), n0, replace=TRUE, CM[1,]),sample(c(0,1,2), n1, replace=TRUE, CM[2,]),sample(c(0,1,2), n2, replace=TRUE, CM[3,]))
    }
  }
  
  
  list(X_0=X_0, X_1=X_1, X_2=X_2, noise_label=noise_label)
}
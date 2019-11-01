pesosKM <- function(time, status){
  
  #time<-db3$Stime
  #status<-db3$event
  M1 <- cbind(time, status) 

  n = nrow(M1)
  
 
  M2 <- matrix(0,nrow=n,ncol=n)
  ord<-order(M1[,2],decreasing=TRUE) 
  
  M2[,2]<-sort(M1[,2],decreasing=TRUE) 
  
  
  M2[,-2]<-M1[ord,-2]  
  
  R=rank(M2[,1],ties.method="first") 
  
  Pkm2=rep(1,n)
  
  Pkm2<-1-M2[,2]/(n-R+1)
  
  count<-outer(R,R,"<") 
  
  
  Pkm2_aux<-matrix(Pkm2,nrow=n,ncol=n,byrow=FALSE) 
 
  
  Pkm2_2<-count*Pkm2_aux
 
  Pkm2_2[Pkm2_2[,]==0]<-1  
  
  Pkm2_cum<-apply(Pkm2_2,2,prod) 
  
  Pkm3<-M2[,2]/(n-R+1) 
  
  Wkm<-rep(0,n) 
  
  Wkm<-Pkm3*Pkm2_cum 
  
  ord2<-order(M2[,1],decreasing=FALSE) 
  
  Wkm<-Wkm[ord2]
  Wkm
  
  ord2<-rank(time,ties.method="first") 
  Wkm<-Wkm[ord2]
  
  table(Wkm==0) 
  
  return(Wkm)
}


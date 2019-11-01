NPPKMW <- function(time, status, bound="none", bw.selec="plug-in", fixed.bw=NULL){
  
  #time=db2$time1
  #status=db2$event1
  #bw.selec=bw.selec
  #fixed.bw=fixed.bw
  #bound=bound
  
  S.pi <- presmooth(times=time, status=status, estimand="S", 
                    bw.selec=bw.selec, bound=bound, fixed.bw = fixed.bw)

  
  x <- S.pi$x.est    
  y <- S.pi$estimate
  
  x1 <- c(0,x)
  y1 <- c(1,y) 
  
  mat <- cbind(x1,y1)
  mat <- unique(mat, MARGIN=1)
  
  dim(mat) 
  length(x1) 
  length(unique(x1))
  length(unique(time))
  
  n <- length(x1) 
  pkw <- rep(0,n-1) 
  
  for(k in 2:n){
    
    pos0 <- which(x1 %in% time[k-1])  #position with repetead times
    
    pos  <- which(mat[,1] %in% time[k-1]) #position when without rep. times
    
    pos1 <- min(pos)
    
    pkw[k-1] <- (mat[pos1-1,2] - mat[pos1,2])/length(pos0)
    
  }
  
  return(pkw)
}
#' Summarizing fits of "pstp" class
#'
#' @description Returns a a data.frame or list containing the estimates of the probabilities.
#' 
#' @param object A fitted pstp object as produced by presmTP.
#' @param state_ini Initial state of the transition. Defaults to state_ini=0.
#' @param times Vector of times; the returned data frame will contain 1 row for each time.
#' @param ... For future methods.


#' @return A data frame or a list containing the estimates of the probability.
#' 
#' @examples
#' res<- presmTP(data = colonIDM, s = 365, method = "uns")
#' summary(res, state_ini=1, times=365*1:5)

#' @author Gustavo Soutinho, Luis Meira-Machado, Pedro Oliveira.

summary.pstp<- function(object, state_ini=0, times = NULL, ...){
  
  #object<-res
  
  #x=object
  
  #times<-time
  
  #object <- x

  if (inherits(object, "pstp")){
    
    
    if (class(object)[1] %in% c("Unsmooth", "Nonparametric", 
                                "Logit", "Logit.gam", "Probit", "Cauchit")) {
      
      
      if (state_ini==0){
        
        if(is.null(times)){
          
          
          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")
          
          res<-object$est0
          print(res)
          
        }else{
          
          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")
          
          resT<-object$est0
          
          #min(resT$t)
          #max(resT$t)
          
          X<-as.data.frame(min(resT$t):max(resT$t))
          colnames(X)<-'t'
          
          dim<-nrow(X)
          
          resM<-merge(x=X, y=resT, by='t',all.x = T)
      
          #head(resM)
          
          for(i in 1:dim){
            #i<-3
            ifelse(is.na(resM$p00[i]),resM$p00[i]<-resM$p00[i-1], resM$p00[i])
            ifelse(is.na(resM$p01[i]),resM$p01[i]<-resM$p01[i-1], resM$p01[i])
            ifelse(is.na(resM$p02[i]),resM$p02[i]<-resM$p02[i-1], resM$p01[i])
          }
          
          res<-resM[resM$t %in% times,]
          print(res)

        }
        
      }else{
        
        if(is.null(times)){
          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")
          res<-object$est1
          print(res)
        }else{
          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")
          #
          #res<-object$est1[object$est1$t %in% times,]
          #print(res)
          #
          resT<-object$est1
          X<-as.data.frame(min(resT$t):max(resT$t))
          colnames(X)<-'t'
          
          dim<-nrow(X)
          
          resM<-merge(x=X, y=resT, by='t',all.x = T)
          
          #head(resM)
          
          for(i in 1:dim){
            #i<-3
            ifelse(is.na(resM$p11[i]),resM$p11[i]<-resM$p11[i-1], resM$p11[i])
            
            ifelse(is.na(resM$p12[i]),resM$p12[i]<-resM$p12[i-1], resM$p12[i])
          
          }

          res<-resM[resM$t %in% times,]
          
          print(res)
          
        }
        
        
      } #fim state_ini
    
      
      
    }else{ #nao ser um dos metodos de estimacao
      
      
      stop("Possible methods are 'uns', 'np', logit, 'logit.gam', 'probit', 'cauchit'.")
    }
    

  }else{
    
       stop("Argument x must be either pstp object.")
    
  }
  
  class(res) <- c("summary.pstp")
  
  return(invisible(res))
  

}




  
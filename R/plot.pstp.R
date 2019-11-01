#' Plot for an object of class "pstp"
#'
#' @description It draws the estimated probabilities.
#' 
#' @param x A fitted pstp object as produced by presmTP.
#' @param state_ini Initial state of the transition. Defaults to state_ini=0
#' @param ... For future methods.

#' @return No value is returned.
#' 
#' @examples
#' res<- presmTP(data = colonIDM, s = 365,method = "uns")
#' plot(res)
#' 
#' @author Gustavo Soutinho, Luis Meira-Machado, Pedro Oliveira.


plot.pstp<- function(x=object,  state_ini=0, ...){
  
  #object<-res
  
  #x=object
  
  ###x<-res #para testar erro!! acrescentei
  
  object <- x

  if (inherits(object, "pstp")){
    
  if (class(object)[1] %in% c("Unsmooth", "Nonparametric", "Logit", "Logit.gam", "Probit", "Cauchit")) {
      
      if(state_ini==0){
        
        times<-object$est0$t
        s<-object$s
        p00<-object$est0$p00
        p01<-object$est0$p01
        p02<-object$est0$p02
        

        x=list(times, cbind(p00, p01, p02))
        
        matplot(x=x[[1]], y=x[[2]], type="l", col=1:3, 
                ylab = paste('pij(', s,',t)', sep=''),lwd = 1)
        
        legend("topright", legend=c("00", "01", "02"), text.col=1:3, cex=1)
       
      }else{
        
        times<-object$est1$t
        s<-object$s
        p11<-object$est1$p11
        p12<-object$est1$p12
        
        x=list(times, cbind(p11, p12))
        
        matplot(x=x[[1]], y=x[[2]], type="l", col=1:2, 
                ylab = paste('pij(', s,',t)', sep=''),lwd = 1)
        
        legend("topright", legend=c("11","12"), text.col=1:2, cex=1)
        
      }

    }else{ #nao ser um dos metodos de estimacao
      
      
      stop("Possible methods are 'uns', 'np', logit, 'logit.gam', 'probit', 'cauchit'.")
    }
    

  }else{
    
       stop("Argument x must be either pstp object.")
    
  }
  
  

} #fim function




  
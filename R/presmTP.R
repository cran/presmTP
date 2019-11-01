#' Methods for estimation of transition probabilities in the illness-death model
#'
#' @description This function is used to obtain unsmoothed and presmoothed estimates of
#' the transition probabilities in the illness-death model.
#' 
#' @param data A numeric value to be squared
#' @param s The first time for obtaining estimates for the transition
#' probabilities.
#' @param method The method used to compute the transition probabilities.
#' Possible options are \code{"uns"}, \code{"np"} \code{"logit"}, \code{"logit.gam"},
#' \code{"probit"} and \code{"cauchit"}. Defaults to \code{"uns"}.
#' @param estimand An optional character string identifying the function to estimate: "S" 
#' for survival function and "H" for cumulative hazard function. Defaults to "S".
#' @param bw.selec An optional (partially matched) character string specifying the method of 
#' bandwidth selection. "fixed" if no bandwidth selection is done, in which case the bandwidth(s)
#' given by the fixed.bw argument is (are) used, "plug-in" for plug-in bandwidth selection and
#' "bootstrap" for bootstrap bandwidth selection. Defaults to "fixed".
#' @param fixed.bw An optional numeric vector with the fixed bandwidth(s) used when the value of
#' the bw.selec argument is "fixed". It must be of length 1 for estimating survival and cumulative
#' hazard functions, and of length 2 for density and hazard functions (in this case, the first element
#' is the presmoothing bandwidth).
#' @param bound An optional numeric vector with the fixed bandwidth(s) used when the value of the bw.selec
#' argument is "fixed". It must be of length 1 for estimating survival and cumulative hazard functions, and
#' of length 2 for density and hazard functions (in this case, the first element is the presmoothing bandwidth).
#' 
#' @return An object of class "pstp" and one of the following classes: \code{"uns"}, 
#' \code{"np"}, \code{"logit"}, \code{"logit.gam"},
#' \code{"probit"} and \code{"cauchit"}. Objects are implemented as a list with elements:
#' \item{est0}{data.frame with estimates of the transition probabilities 0->0, 0->1 and 0->2.}
#' \item{est1}{data.frame with estimates of the transition probabilities 1->1 and 1->2.}
#' \item{s}{The first time for obtaining estimates for the transition probabilities.}
#' \item{callp}{The expression of the estimated probability.}
#' \item{call}{A call object.}
#' 
#' @examples
#' 
#' #Unsmoothed
#' res1<- presmTP(data = colonIDM, s = 365,method = "uns" )
#' res1$est0$t
#' res1$est0$p02
#' res1$est1$t

#' summary(res1, state_ini=1, time=365*1:5)
#' plot(res1)

#' res1$call
#' class(res1)

#' #Nonparametric
#' res2<- presmTP(data = colonIDM, s = 365,method = "np" )
#' res3<- presmTP(data = colonIDM, s = 365,method = "np", estimand="S")
#' res4<- presmTP(data = colonIDM, s = 365,method = "np", estimand="H")
#' res5<- presmTP(data = colonIDM, s = 365,method = "np", 
#'                bw.selec="fixed", fixed.bw=30)

#' #Presmoothed - Logit
#' res6<- presmTP(data = colonIDM, s = 365,method = "logit" )
#' summary(res6, state_ini=1, time=365*1:5)

#' #Presmoothed - Logit GAM
#' res7<- presmTP(data = colonIDM, s = 365,method = "logit.gam" )
#' 
#' 
#' 
#' @references
#' 
#' Aalen O. O., Johansen S. (1978) An Empirical Transition Matrix for
#' Nonhomogeneous Markov Chains Based on Censored Observations. Scandinavian
#' Journal of Statistics 5(3), 141--150.
#' 
#' Cao, R., Lopez-de-Ullibarri, I., Janssen, P. and Veraverbeke, N. (2005). Presmoothed Kaplan-Meier 
#' and Nelson-Aalen estimators, Journal of Nonparametric Statistics, 17, 31-56.
#'
#' Meira-Machado L. F., de Una-Alvarez J. and Cadarso-Suarez C. (2006).
#' Nonparametric estimation of transition probabilities in a non-Markov
#' illness-death model. Lifetime Data Anal 12(3), 325--344.
#' 
#' Lopez-de-Ullibarri, I and Jacome, M. A. (2013). survPresmooth: An R Package for Presmoothed
#' Estimation in Survival Analysis, Journal of Statistical Software, 54(11), 1-26. 
#' URL: http://www.jstatsoft.org/v54/i11/.

#' de Una-Alvarez J. and Meira-Machado L. (2015). Nonparametric estimation
#' of transition probabilities in a non-Markov illness-death model:
#' a comparative study. Biometrics 71, 364--375.
#
#' 
#' Meira-Machado, L. (2016). Smoothed landmark estimators of the transition probabilities, 
#' SORT-Statistics and Operations Research Transactions, 40, 375-398.
#' 


#'@author Gustavo Soutinho, Luis Meira-Machado, Pedro Oliveira.



#' @importFrom "survPresmooth" presmooth
#' @importFrom "mgcv" gam
#' @importFrom "graphics" legend matplot plot
#' @importFrom "stats" binomial glm predict
#' @importFrom "utils" head

#' @export presmTP
#' @export summary.pstp
#' @export plot.pstp
#' 
#' @S3method summary pstp
#' @S3method plot pstp
#' 
#' 
presmTP<-function(data, s, method = "uns", estimand="S", 
                  bw.selec="plug-in", fixed.bw=NULL, bound="none"){
  
  if (missing(s))
    stop("argument 's' is missing, with no default")
  
  if (!(method %in% c("uns", "np", "logit", "logit.gam", "probit", "cauchit"))) {
    stop("Possible methods are 'uns', 'np', logit, 'logit.gam', 'probit', 'cauchit'.")
  }
  
    obj_data <- data
    
    p <- which(obj_data$time1 > s) 
    
    p 
    
    db2<- obj_data [p,] 
    
    times <- c(db2$Stime, db2$time1, s) 
    
    times <- unique(sort(times)) 
    
    n <- length(times); n #788
    
    if(method=='uns'){
      
      p11.unsm <- rep(0,n) 
      
      resZunsmooth <- pesosKM(db2$time1, db2$event1)
      table(resZunsmooth==0) 
      
      
      for(k in 1:n){
        
        p <- which(db2$time1 <= times[k])
        
        p11.unsm[k] <- 1-sum(resZunsmooth[p]) #Sobrevivencia=1-sum(pesos) (z)
        
      }
      
      
      p13.unsm <- rep(0,n)
      
      resTunsm <- pesosKM(db2$Stime, db2$event) 
      
      
      for(k in 1:n){
        
        p <- which(db2$Stime <= times[k])
        
        p13.unsm[k] <- sum(resTunsm[p]) 
        
      }
      
      
      p12.unsm <- rep(0,n)
      
      for(k in 1:n){
        
        p1 <- which(db2$time1 <= times[k])
        p2 <- which(db2$Stime <= times[k])
        
        p12.unsm[k] <- sum(resZunsmooth[p1]) - sum(resTunsm[p2])
        
      }
      
      
      p <- which(obj_data$time1 <= s & obj_data$Stime > s)
      
      length(p) 
      
      db3<-obj_data[p,] #z<=s ^  T>s
      
      times2 <- c(db3$Stime, db3$time1,s)
      times2 <-times2[times2 >= s]
      times2 <- unique(sort(times2)) 
      m <- length(times2); m  
      
      p22.unsm <- rep(0,m)
      
      res22unsm <- pesosKM(db3$Stime, db3$event)
      
      #length(res22unsm) #152
      
      for(k in 1:m){
        p <- which(db3$Stime <= times2[k])
        p22.unsm[k] <- 1 - sum(res22unsm[p])
        
        
      }
      
      
    p23.unsm <- rep(0,m)
      
      
      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        p23.unsm[k] <-sum(res22unsm[p])
        
        
      }
      
      
      #
      
      resu0 <- data.frame(cbind(times, p11.unsm, p12.unsm, p13.unsm))  
      
      names(resu0) <- c("t", "p00", "p01", "p02")
      
      head(resu0)
      
      resu1 <- data.frame(cbind(times2, p22.unsm, p23.unsm))  
      
      names(resu1) <- c("t", "p11", "p12")
      
      head(resu1)
      
      res <- list(est0 = resu0,est1= resu1,s = s)  
      
      class(res) <- c("Unsmooth", "pstp")
      
      
      
    }
    
    
    if(method=='np' & estimand=='S'){
      
      
      p11.np <- rep(0,n)  
      
      
      #Nonparametric presmoothing:
      resZnp <- NPPKMW(db2$time1, db2$event1,bw.selec=bw.selec, fixed.bw=fixed.bw, bound=bound) 
      #table(resZnp==0) 
      
      
      for(k in 1:n){
        
        p <- which(db2$time1 <= times[k])
        p11.np[k] <- 1 - sum(resZnp[p])
        
      }
      
      
      p13.np <- rep(0,n)     #nonparametric presmoothing
      
      resTnp <- NPPKMW(db2$Stime, db2$event,bw.selec=bw.selec, fixed.bw=fixed.bw, bound=bound) 
      
      for(k in 1:n){
        
        p <- which(db2$Stime <= times[k])
        
        p13.np[k] <- sum(resTnp[p]) 
        
      }
      
      p12.np <- rep(0,n)   
      
      
      for(k in 1:n){
        
        p1 <- which(db2$time1 <= times[k])
        p2 <- which(db2$Stime <= times[k])
        
        p12.np[k] <- sum(resZnp[p1]) - sum(resTnp[p2])
        
      }
      
      p <- which(obj_data$time1 <= s & obj_data$Stime > s)
      
      length(p) 
      
      db3<-obj_data[p,] #z<=s ^  T>s
      
      times2 <- c(db3$Stime, db3$time1,s)
      times2 <-times2[times2 >= s]
      times2 <- unique(sort(times2)) 
      m <- length(times2); m  
      
      p22.np <- rep(0,m)     
      
       
      res22np <-suppressWarnings(tryCatch({NPPKMW(db3$Stime, db3$event,bw.selec=bw.selec, fixed.bw=fixed.bw, bound=bound)}, 
                          error=function(j){('erro22')}))
      
      if(unique(res22np=='erro22')){
        
        warning("Presmoothed estimates of transition 1->1 computed with the fixed bandwidth 0.01. 
        Default method 'plug-in' unable to be applied.")
        
        res22np<-NPPKMW(db3$Stime, db3$event, bw.selec='fixed', fixed.bw=0.01, bound=bound)

        for(k in 1:m){
          
          p <- which(db3$Stime <= times2[k])
          
          p22.np[k] <- 1 - sum(res22np[p])
          
        }
        
      }else{
        
        for(k in 1:m){
          
          p <- which(db3$Stime <= times2[k])
          
          p22.np[k] <- 1 - sum(res22np[p])

        }
        
      }
      
      
      p23.np <- rep(0,m)   
      
      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        
        p23.np[k] <-sum(res22np[p])
        
      }
      
      
      #
      
      resu0 <- data.frame(cbind(times, p11.np, p12.np, p13.np))  
      
      names(resu0) <- c("t", "p00", "p01", "p02")
      
      head(resu0)
      
      resu1 <- data.frame(cbind(times2, p22.np, p23.np))  
      
      names(resu1) <- c("t", "p11", "p12")
      
      head(resu1)
      
      res<- list(est0 = resu0,est1 = resu1,s = s)  
      
      class(res) <- c("Nonparametric", "pstp")
      
      #return(invisible(result))
      
    }
    
    if(method=='np' & estimand=='H'){
      
      
      p11.fh <- rep(0,n)  
      
      #Nonparametric presmoothing:
      resZfh <- FH(db2$time1, db2$event1,bw.selec=bw.selec, fixed.bw=fixed.bw, bound=bound) 
      table(resZfh==0) 
      
      
      for(k in 1:n){
        
        p <- which(db2$time1 <= times[k])
        
        p11.fh[k] <- 1 - sum(resZfh[p])
        
        
      }
      
      
      p13.fh <- rep(0,n)     #nonparametric presmoothing
      
      resTfh <- FH(db2$Stime, db2$event,bw.selec=bw.selec, fixed.bw=fixed.bw, bound=bound) 
      
      for(k in 1:n){
        
        p <- which(db2$Stime <= times[k])
        
        p13.fh[k] <- sum(resTfh[p]) 
        
      }
      
      p12.fh <- rep(0,n)   
      
      
      for(k in 1:n){
        
        p1 <- which(db2$time1 <= times[k])
        p2 <- which(db2$Stime <= times[k])
        
        p12.fh[k] <- sum(resZfh[p1]) - sum(resTfh[p2])
        
      }
      
      p <- which(obj_data$time1 <= s & obj_data$Stime > s)
      
      length(p) 
      
      db3<-obj_data[p,] #z<=s ^  T>s
      
      times2 <- c(db3$Stime, db3$time1,s)
      times2 <-times2[times2 >= s]
      times2 <- unique(sort(times2)) 
      m <- length(times2); m  

      p22.fh <- rep(0,m)     
      
      res22fh <-   suppressWarnings(FH(db3$Stime, db3$event,bw.selec=bw.selec, fixed.bw=fixed.bw, bound=bound))
      
      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        
        p22.fh[k] <- 1 - sum(res22fh[p])
 
      }
      p23.fh <- rep(0,m)   
      
      
      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        
        p23.fh[k] <-sum(res22fh[p])
        
        
      }

      #
      resu0 <- data.frame(cbind(times, p11.fh, p12.fh, p13.fh))  
      
      names(resu0) <- c("t", "p00", "p01", "p02")
      
      head(resu0)
      
      resu1 <- data.frame(cbind(times2, p22.fh, p23.fh))  
      
      names(resu1) <- c("t", "p11", "p12")
      
      head(resu1)
      
      res<- list(est0 = resu0,est1 = resu1,s = s)  
      
      class(res) <- c("Nonparametric", "pstp")
      
      #return(invisible(result))
      
    }
    
    
    if(method=='logit'){
      
      p11.logit <- rep(0,n)
      
      
      #Deltas presmoothing
      
      fit.logit0 <- glm(event1 ~ time1, data = db2, family=binomial(link="logit"))
      prob.logit0 <- predict(fit.logit0, type="response")
      
      resZlogit <- pesosKM(db2$time1, prob.logit0) #call function
      table(resZlogit ==0) 
      
      
      
      for(k in 1:n){
        
        p <- which(db2$time1 <= times[k])
        
        
        p11.logit[k] <- 1 - sum(resZlogit[p]) 
        
      }
      
      
      p13.logit <- rep(0,n)
      
      fit.logit <- glm(event ~ Stime, data = db2, family=binomial(link="logit"))
      prob.logit <-predict(fit.logit,type="response")
      resTlogit <- pesosKM(db2$Stime, prob.logit)
      
      for(k in 1:n){
        
        p <- which(db2$Stime <= times[k])
        
        
        p13.logit[k] <- sum(resTlogit[p])
        
      }
      
      p12.logit <- rep(0,n)
      
      for(k in 1:n){
        
        p1 <- which(db2$time1 <= times[k])
        p2 <- which(db2$Stime <= times[k])
        
        p12.logit[k] <- sum(resZlogit[p1]) - sum(resTlogit[p2])
        
      }
      
      p <- which(obj_data$time1 <= s & obj_data$Stime > s)
      
      length(p) 
      
      db3<-obj_data[p,] #z<=s ^  T>s
      
      times2 <- c(db3$Stime, db3$time1,s)
      times2 <-times2[times2 >= s]
      times2 <- unique(sort(times2)) 
      m <- length(times2); m  
      
      
      p22.logit <- rep(0,m)
      
      fit.logit <- glm(event ~ Stime, data = db3, family=binomial(link="logit"))
      prob.logit <-predict(fit.logit,type="response")
      
      res22logit <- pesosKM(db3$Stime, prob.logit)
      
      for(k in 1:m){
        p <- which(db3$Stime <= times2[k])
        
        p22.logit[k] <- 1 - sum(res22logit[p])
        
        
      }
      p23.logit <- rep(0,m)

      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        
        p23.logit[k] <- sum(res22logit[p])
        
        
      }
      
      
      #
      
      resu0 <- data.frame(cbind(times, p11.logit, p12.logit, p13.logit))  
      
      names(resu0) <- c("t", "p00", "p01", "p02")
      
      head(resu0)
      
      resu1 <- data.frame(cbind(times2, p22.logit, p23.logit))  
      
      names(resu1) <- c("t", "p11", "p12")
      
      head(resu1)
      
      res<- list(est0 = resu0,est1 = resu1,s = s)  
      
      class(res) <- c("Logit", "pstp")
      
      #return(invisible(result))
      
    }
    
    if(method=='logit.gam'){
      
      p11.logit3<- rep(0,n) #using gam
      
      fit.logitgam <- gam(event1 ~ s(time1), data = db2, family=binomial(link="logit")) #CALL "mgcv"
      prob.logitgam <-predict(fit.logitgam,type="response")
      
      
      resZlogit3 <- pesosKM(db2$time1, prob.logitgam)
      
      
      
      for(k in 1:n){
        
        p <- which(db2$time1 <= times[k])
        
        p11.logit3[k] <- 1 - sum(resZlogit3[p]) 
        
        
        
      }
      
      
      p13.logit3 <- rep(0,n) #using gam
      
      fit.logitgam <- gam(event ~ s(Stime), data = db2, family=binomial(link="logit"))
      prob.logitgam <-predict(fit.logitgam,type="response")
      
      
      resTlogit3 <- pesosKM(db2$Stime, prob.logitgam)
      
      for(k in 1:n){
        
        p <- which(db2$Stime <= times[k])
        
        p13.logit3[k] <- sum(resTlogit3[p])
        
      }
      
      p12.logit3 <- rep(0,n) #using gam
      
      
      for(k in 1:n){
        
        p1 <- which(db2$time1 <= times[k])
        p2 <- which(db2$Stime <= times[k])
        
        p12.logit3[k] <- sum(resZlogit3[p1]) - sum(resTlogit3[p2])
        
        
      }
      
      p <- which(obj_data$time1 <= s & obj_data$Stime > s)
      
      length(p) 
      
      db3<-obj_data[p,] #z<=s ^  T>s
      
      times2 <- c(db3$Stime, db3$time1,s)
      times2 <-times2[times2 >= s]
      times2 <- unique(sort(times2)) 
      m <- length(times2); m  
      
      p22.logit3 <- rep(0,m) #using gam
      
      
      fit.logitgam <- gam(event ~ s(Stime), data = db3, family=binomial(link="logit"))
      prob.logitgam <-predict(fit.logitgam,type="response")
      
      res22logit3 <- pesosKM(db3$Stime, prob.logitgam)
      
      for(k in 1:m){
        p <- which(db3$Stime <= times2[k])
        
        p22.logit3[k] <- 1 - sum(res22logit3[p])
        
        
      }
      p23.logit3 <- rep(0,m) #using gam
      

      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        
        p23.logit3[k] <- sum(res22logit3[p])
        
        
      }
      
      #
      resu0 <- data.frame(cbind(times, p11.logit3, p12.logit3, p13.logit3))  
      
      names(resu0) <- c("t", "p00", "p01", "p02")
      
      head(resu0)
      
      resu1 <- data.frame(cbind(times2, p22.logit3, p23.logit3))  
      
      names(resu1) <- c("t", "p11", "p12")
      
      head(resu1)
      
      res<- list(est0 = resu0,est1 = resu1,s = s)  
      
      class(res) <- c("Logit.gam", "pstp")
      
      #return(invisible(result))
      
    }
    
    
    if(method=='probit'){
      
      p11.probit <- rep(0,n)
      
      fit.probit0 <- glm(event1 ~ time1, data = db2, family=binomial(link="probit"))
      prob.probit0 <- predict(fit.probit0, type="response")
      
      resZprobit <- pesosKM(db2$time1, prob.probit0)
      
      
      
      for(k in 1:n){
        
        p <- which(db2$time1 <= times[k])
        
        p11.probit[k] <- 1 - sum(resZprobit[p])
        
        
      }
      
      p13.probit <- rep(0,n)
      
      
      fit.probit <- glm(event ~ Stime, data = db2, family=binomial(link="probit"))
      prob.probit <-predict(fit.probit,type="response")
      
      resTprobit <- pesosKM(db2$Stime, prob.probit)
      
      for(k in 1:n){
        
        p <- which(db2$Stime <= times2[k])
        
        p13.probit[k] <- sum(resTprobit[p])
        
        
      }
      
      
      p12.probit <- rep(0,n)
      
      for(k in 1:n){
        
        p1 <- which(db2$time1 <= times[k])
        p2 <- which(db2$Stime <= times[k])
        
        p12.probit[k] <- sum(resZprobit[p1]) - sum(resTprobit[p2])
        
        
      }
      
      p <- which(obj_data$time1 <= s & obj_data$Stime > s)
      
      length(p) 
      
      db3<-obj_data[p,] #z<=s ^  T>s
      
      times2 <- c(db3$Stime, db3$time1,s)
      times2 <-times2[times2 >= s]
      times2 <- unique(sort(times2)) 
      m <- length(times2); m  
      
      p22.probit <- rep(0,m)
      
      fit.probit <- glm(event ~ Stime, data = db3, family=binomial(link="probit"))
      prob.probit <-predict(fit.probit,type="response")
      res22probit <- pesosKM(db3$Stime, prob.probit)
      
      for(k in 1:m){
        p <- which(db3$Stime <= times2[k])
        
        p22.probit[k] <- 1 - sum(res22probit[p])
  
      }
      p23.probit <- rep(0,m)
      

      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        
        p23.probit[k] <-sum(res22probit[p])
        
        
      }
      
      
      #
      
      resu0 <- data.frame(cbind(times, p11.probit, p12.probit, p13.probit))  
      
      names(resu0) <- c("t", "p00", "p01", "p02")
      
      head(resu0)
      
      resu1 <- data.frame(cbind(times2, p22.probit, p23.probit))  
      
      names(resu1) <- c("t", "p11", "p12")
      
      head(resu1)
      
      res<- list(est0 = resu0,est1 = resu1,s = s)  
      
      class(res) <- c("Probit", "pstp")
      
      #return(invisible(result))
      
    }
    
    if(method=='cauchit'){
      
      p11.cauchit <- rep(0,n)
      
      fit.cauchit0 <- glm(event1 ~ time1, data = db2, family=binomial(link="cauchit"))
      prob.cauchit0 <-predict(fit.cauchit0, type="response")
      
      resZcauchit <- pesosKM(db2$time1, prob.cauchit0)
      
      
      
      for(k in 1:n){
        
        p <- which(db2$time1 <= times[k])
        
        p11.cauchit[k] <- 1 - sum(resZcauchit[p]) 
        
      }
      
      p13.cauchit <- rep(0,n)
      
      
      fit.cauchit <- glm(event ~ Stime, data = db2, family=binomial(link="cauchit"))
      prob.cauchit <-predict(fit.cauchit,type="response")
      
      
      resTcauchit <- pesosKM(db2$Stime, prob.cauchit)
      
      for(k in 1:n){
        
        p <- which(db2$Stime <= times[k])
        
        p13.cauchit[k] <- sum(resTcauchit[p])
        
      }
      

      p12.cauchit <- rep(0,n)
      
      
      for(k in 1:n){
        
        p1 <- which(db2$time1 <= times[k])
        p2 <- which(db2$Stime <= times[k])
        
        p12.cauchit[k] <- sum(resZcauchit[p1]) - sum(resTcauchit[p2])
        
      }
      
      p <- which(obj_data$time1 <= s & obj_data$Stime > s)
      
      length(p) 
      
      db3<-obj_data[p,] #z<=s ^  T>s
      
      times2 <- c(db3$Stime, db3$time1,s)
      times2 <-times2[times2 >= s]
      times2 <- unique(sort(times2)) 
      m <- length(times2); m  #141
      
      p22.cauchit <- rep(0,m)
      
      fit.cauchit <- glm(event ~ Stime, data = db3, family=binomial(link="cauchit"))
      prob.cauchit <-predict(fit.cauchit,type="response")
      res22cauchit <- pesosKM(db3$Stime, prob.cauchit)
      
      for(k in 1:m){
        p <- which(db3$Stime <= times2[k])
        
        p22.cauchit[k] <- 1 - sum(res22cauchit[p])
        
      }
      p23.cauchit <- rep(0,m)
      
      
      
      
      for(k in 1:m){
        
        p <- which(db3$Stime <= times2[k])
        
        p23.cauchit[k] <-sum(res22cauchit[p])
        
      }
      #
      resu0 <- data.frame(cbind(times, p11.cauchit, p12.cauchit, p13.cauchit))  
      
      names(resu0) <- c("t", "p00", "p01", "p02")
      
      head(resu0)
      
      resu1 <- data.frame(cbind(times2, p22.cauchit, p23.cauchit))  
      
      names(resu1) <- c("t", "p11", "p12")
      
      head(resu1)
      
      res<- list(est0 = resu0,est1 = resu1,s = s)  
      
      class(res) <- c("Cauchit", "pstp")
      
      #return(invisible(result))
      
    }
    suppressWarnings(res)
    callp <- paste("pij(s=", s, ",t)", sep = "")
    res$callp <- callp
    res$call <- match.call()
   
    return(invisible(res))


}



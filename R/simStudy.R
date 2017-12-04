#' @export
#' @importFrom utils head tail 
#' @importFrom stats rbinom rgamma rpois na.omit aggregate 

simStudy <- function(n,numsim,beta,nu,mu,inf=FALSE,nobs){
  
  simResults <- list()
  estimates <- matrix(rep(NA,4*numsim),ncol=4)
  
  for(i in 1:numsim){
    set.seed(i)
    
    if(inf==FALSE){
      simdata <- simRtc(n=n,beta=beta,nu=nu,mu=mu,inf=FALSE,nobs=nobs)
    }
    
    if(inf==TRUE){
      simdata <- simRtc(n=n,beta=beta,nu=nu,mu=mu,inf=TRUE,nobs=c(nobs[1],nobs[2]))
    }
    
    print(paste('Generated data for dataset ',i))
    
    gee.fit <- suppressMessages(gee(min ~ x1 ,data = simdata, id = ID, family = "poisson", corstr = "independence",silent=T))
    print(paste('GEE fit for dataset ',i))
    
    fitted.gee <- c(summary(gee.fit)$coefficients[2,1],
                    summary(gee.fit)$coefficients[2,4])
    names(fitted.gee) <- c('estimate','se')
    print(fitted.gee)
    
    formula <- aee(ID, time, min) ~ x1 #+ x2
    fitted2 <- aeexfit(formula=formula, data=simdata,se="Sandwich",weight=NULL) 
    fitted.aeex <- c(fitted2$beta,fitted2$betaSE)
    print(paste('Semiparametric model fit for dataset ',i))
    names(fitted.aeex) <- c('estimate','se')
    print(fitted.aeex)
    
    estimates[i,] <- c(fitted.gee,fitted.aeex)
    
    colnames(estimates) <- c('gee.estimate','gee.se','aeex.estimate','aeex.se')
  }
  
  simResults$estimates <- estimates
  colnames(simResults$estimates) <- c('gee.estimate','gee.se','aeex.estimate','aeex.se')
  
  bias <- rep(NA,2)
  bias[1] <- mean(estimates[,1],na.rm=TRUE) - (beta)
  bias[2] <- mean(estimates[,3],na.rm=TRUE) - (beta)
  
  simResults$bias <- bias
  names(simResults$bias) <- c('gee.bias','aeex.bias')
  
  covered.gee <- ifelse(estimates[,1]-1.96*estimates[,2] <= beta &
                          estimates[,1]+1.96*estimates[,2] >= beta, 1, 0)
  
  covered.aee <- ifelse(estimates[,3]-1.96*estimates[,4] <= beta &
                          estimates[,3]+1.96*estimates[,4] >= beta, 1, 0)
  
  simResults$coverage <- c(mean(covered.gee,na.rm=TRUE),mean(covered.aee,na.rm=TRUE))
  names(simResults$coverage) <- c('gee.coverage','aeex.coverage')
  
  simResults
}

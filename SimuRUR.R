# ---------------------------------------------------------------------------
# FUNÇÃO PARA SIMULAÇÃO DE SÉRIES TEMPORAIS USANDO O MODELO RUR ARMA
# ---------------------------------------------------------------------------
simu.RUR <- function(n,phi=0.2,theta=0.4, alpha=1, tau=0.5, freq=12,
                             link="logit") # phi = Autoregressivo             theta = Médias moveis           tau = quantil
{
  source("Funcoes.R")
  
  
  if(any(is.na(phi)==F))
  {
    ar <- 1:length(phi)
  }else{
    ar<-0
    phi<-0
  }
  
  if(any(is.na(theta)==F))
  {
    ma <- 1:length(theta)
  }else{
    ma<-0
    theta<-0
  }
  
  linktemp <- substitute(link)
  if (!is.character(linktemp))
  {
    linktemp <- deparse(linktemp)
    if (linktemp == "link")
      linktemp <- eval(link)
  }
  if (any(linktemp == c("logit", "probit", "cloglog")))
  {  
    stats <- make.link(linktemp)
  }else{
    stop(paste(linktemp, "link not available, available links are \"logit\" 
               and \"cloglog\""))
  } 
  
  link <- structure(list(link = linktemp, 
                         linkfun = stats$linkfun,
                         linkinv = stats$linkinv ))
  
  linkfun <- link$linkfun
  linkinv <- link$linkinv
  
  {
    p <- max(ar)
    q <- max(ma)
    m <- 2*max(p,q)
    
    ynew <-rep(alpha,(n+m))
    mu <- linkinv(ynew)
    
    error<-rep(0,n+m) 
    eta<- y <- rep(NA,n+m)
    
    for(i in (m+1):(n+m))
    {
      eta[i]  <- alpha + as.numeric(phi%*%ynew[i-ar]) + as.numeric(theta%*%error[i-ma])
      mu[i]   <- linkinv(eta[i])
      #u <- runif(1)
      y[i]    <- rRUR(1,mu[i],tau)
      # if(y[i]==1){
      #   u <- runif(1)
      #   y[i]    <- qLB(u,mu[i])
      ynew[i] <- linkfun(y[i])
      error[i]<- ynew[i]-eta[i]
    }
    
    
    return(ts(y[(m+1):(n+m)],frequency=freq) )
  } 
}


# plot(simu.RUR(100))
# y<-simu.RUR(10000)
# plot(simu.RUR(1000))

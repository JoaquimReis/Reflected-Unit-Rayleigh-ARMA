# density function Função densidade
dRUR<-function(y, mu = 0.7, tau = 0.5, log = FALSE)
{
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(y <= 0) | any(y >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  fy1 <- ((2*log(1-tau))/(y-1))*(-log(1-mu))^(-2)*(-log(1-y))^(2-1)*(1-tau)^((-log(1-y)/-log(1-mu))^2)
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}

integrate(dRUR,0,1)
#------------------------------------------------------------------------------------------ #ok
# cumulative distribution function Acumulada
pRUR<-function(q, mu = 0.7, lower.tail = TRUE, log.p = FALSE){
  if (any(mu <= 0) | any(mu >= 1)) stop(paste("mu must be between 0 and 1", "\n", ""))
  if (any(q <= 0) | any(q >= 1)) stop(paste("y must be between 0 and 1", "\n", ""))
  cdf1 <- 1-(1-0.5)^(((-log(1-q))/(-log(1-mu)))^2)
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}

pRUR(.5)
integrate(dRUR,0,.5)
#------------------------------------------------------------------------------------------ #ok
# quantile function Quantil
qRUR<-function(u, mu)
{
  q <- 1-(1-mu)^(-log(1-u)/-log(1-0.5))^(1/2)
  q
}

u=pRUR(.5)
qRUR(u,mu=.7)

# inversion method for random generation Gerador de amostras aleatórias
rRUR<-function(n,mu,a=0,b=1)
{
  u<- runif(n)
  tau<-.5
  y <- qRUR(u,mu = mu)
  y <- pmin(pmax(y, 1e-6), 1 - 1e-6)
  y
}

rRUR(500, 0.4)

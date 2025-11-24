# Exemplo de cálculo para o vetor escore da Kumaraswamy
# parametros q e phi, de quantil e forma respectivamente
# Renata Rojas Guerra (renata.r.guerra@ufsm.br), january/2023
# Reflected Unit Rayleigh

# função densidade
druw<-function(y,mu,tau=.5)
{
  d<-((2*log(1-tau))/(y-1))*(-log(1-mu))^(-2)*
    (-log(1-y))*(1-tau)^((log(1-y)/log(1-mu))^2)
  d
}



#calculando a log-verossimilhança analítica
loglik<-expression(log(((2*log(1-tau))/(y-1))*(-log(1-mu))^(-2)*
                         (-log(1-y))*(1-tau)^((log(1-y)/log(1-mu))^2)))



# definindo valores para testar
mu=.70;tau=.50;y<-.25

((2*log(1-tau))/(y-1))*(-log(1-mu))^(-2)*
  (-log(1-y))*(1-tau)^((log(1-y)/log(1-mu))^2)


# comparando a loglik analitica e numerica
log(druw(y,mu)) # numerica
eval(loglik) # analitica

# calculando a derivada da loglik com respeito a q
D(loglik,"mu") # analitica do R


#avaliando nos valores definidos para testar
eval(D(loglik,"mu"))



# analitica simplificada
(2/((1-mu)*log(1-mu)))*(1+log(1-tau)*(log(1-y)/log(1-mu))^(2))

# NÃO USAR, TESTES analitica simplificada NÃO USAR, TESTES
1/2-log(-log(1-mu))+log(-log(1-y))+log(1-mu)^-2*log(1-tau)*log(1-y)^2*log(log(1-y)/log(1-mu))

1/2-log(-log(1-mu))+log(-log(1-y))+(log(1-y)/log(1-mu))^2*log(1-tau)*log(log(1-y)/log(1-mu))



#analítica dmu 
dmu <- expression((2/((1-mu)*log(1-mu)))*(1+log(1-tau)*(log(1-y)/log(1-mu))^(2)))


#númerica dmu2
D(dmu, "mu")

#avaliando nos valores definidos para testar
eval(D(dmu,"mu"))


#analítica dmu2
(2*(1 + log(1 - mu)))/((-1 + mu)^2*log(1 - mu)^2) + (2*(1 + 2 + log(1 - mu))*log(1 - tau)*(log(1 - y)/log(1 - mu))^2)/((-1 + mu)^2 *log(1 - mu)^2)

#numérica do R 
2 * ((1 - mu) * (1/(1 - mu)) + 
           log(1 - mu))/((1 - mu) * log(1 -  mu))^2 * (1 + log(1 - tau) * (log(1 - y)/log(1 - mu))^(2)) + 
  (2/((1 - mu) * log(1 - mu))) * (log(1 - tau) * ((log(1 -  y)/log(1 - mu))^((2) - 1) * ((2) * (log(1 - y) * 
                                                                                                              (1/(1 - mu))/log(1 - mu)^2))))

#Esperanças numéricas

qq1ruw<-function(y, mu = 0.7, tau = 0.5, log = FALSE)
{
  d<-((log(1-y)/log(1-mu))^2)*((2*log(1-tau))/(y-1))*(-log(1-mu))^(-2)*
    (-log(1-y))*(1-tau)^((log(1-y)/log(1-mu))^2)
  d
}


integrate(qq1ruw, 0, 1)

qq2ruw<-function(y, mu = 0.7, tau = 0.5, log = FALSE)
{
  d<-log((log(1-y)/log(1-mu)))*(log(1-y)/log(1-mu))^2*((2*log(1-tau))/(y-1))*(-log(1-mu))^(-2)*
    (-log(1-y))*(1-tau)^((log(1-y)/log(1-mu))^2)
  d
}


integrate(qq2ruw, 0, 1)

qq3ruw<-function(y, mu = 0.7, tau = 0.5, log = FALSE)
{
  d<-(log((log(1-y)/log(1-mu)))^2)*(log(1-y)/log(1-mu))^2*((2*log(1-tau))/(y-1))*(-log(1-mu))^(-2)*
    (-log(1-y))*(1-tau)^((log(1-y)/log(1-mu))^2)
  d
}


integrate(qq3ruw, 0, 1)



#Esperanças analíticas
EulerGamma <- 0.5772156649

#qq1
-1/log(1-tau)

#qq3
(EulerGamma+log(-log(1-tau))-1)/(2*log(1-tau))

#qq2
numerador <- 6 * log(-log(1-tau)) * (log(-log(1-tau)) + 2 * EulerGamma - 2) + 6 * EulerGamma * (EulerGamma - 2) + pi^2
denominador <- -(6 * log(1-tau) * 2^2)
resultado <- numerador / denominador
resultado


(6 * log(-log(1-tau)) * (log(-log(1-tau)) + 2 * EulerGamma - 2) + 6 * EulerGamma * (EulerGamma - 2) + pi^2) / (-(6 * log(1-tau) * sigma^2))


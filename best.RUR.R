best.RUR <- function(serie, sf, h=6, pmax=6, qmax=6, nbest=10,
                     tau=0.5, link="logit", X=NA, X_hat=NA)
{
  source("Funcoes.R")
  source("FitRUR.R")
  
  y <- ts(serie, start=c(sf[1], sf[2]), frequency=sf[3])
  
  # It initializes the AIC criteria
  fit    <- RURarma.fit(y, ma=1, diag=0, link=link)
  aicmin <- fit$aic
  
  print(aicmin)
  
  best_aic <- rep(Inf, nbest)                    # It saves the nbest smallest AICs
  melhores <- matrix(rep(0,(nbest*3)), ncol=3)   # It saves the order of the nbest best models
  colnames(melhores) <- c("p","q","AIC")
  
  tot <- 0
  bug <- 0
  
  for(p in 0:pmax)
  {
    for(q in 0:qmax)
    {
      if(p==0) ar1 <- NA else ar1 <- 1:p
      if(q==0) ma1 <- NA else ma1 <- 1:q
      
      # pula (0,0): RURarma.fit não aceita ar=NA e ma=NA simultaneamente
      if(sum(is.na(c(ar1,ma1))) < 2)
      {
        print(c(p,q), quote=F)
        fitRUR <- RURarma.fit(y, ar=ar1, ma=ma1, tau=tau, link=link,
                              X=X, X_hat=X_hat, h=h)
        tot <- tot + 1
        
        if(fitRUR$conv != 0)
        {
          print(c("NO CONVERGENCE  ",p,q), quote=F)
          bug <- bug + 1
          next
        }
        
        if(aicmin > fitRUR$aic)   # best model according to AIC
        {
          aicmin         <- fitRUR$aic
          best_model_aic <- fitRUR$model
          best_p         <- p          # <— guarda ordem do melhor modelo
          best_q         <- q
          print("###########################################")
          print(aicmin)
        }
        
        if(fitRUR$aic < max(best_aic))
        {
          maximo            <- order(best_aic)[nbest]
          best_aic[maximo]  <- fitRUR$aic
          melhores[maximo,] <- c(p, q, fitRUR$aic)
        }
      }
    }
  }
  
  print(" ", quote=F)
  print("SELECTED MODEL FROM AIC", quote=F)
  print(best_model_aic, quote=F)
  print(" ", quote=F)
  print(c("Total of tested models =", tot), quote=F)
  print(c("Total of errors in the estimation =", bug), quote=F)
  print(" ", quote=F)
  print("THE BEST MODELS", quote=F)
  print(melhores, quote=F)
  
  # retorna a melhor ordem para uso no script principal
  # análogo ao que orbxiiarma <- c(1,2) faz manualmente no original
  invisible(list(best_order = c(best_p, best_q),
                 best_aic   = aicmin,
                 best_model = best_model_aic,
                 all_best   = melhores))
}

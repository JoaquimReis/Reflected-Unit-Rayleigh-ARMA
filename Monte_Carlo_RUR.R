
#SIMULACOA DE MONTE CARLO ORIGINAL

# simu - ARMA(1,1) - apenas medias moveis

rm(list = ls())

source("simuRUR.R")
source("FitRUR.R")

alpha = .7
phi =  0.4  #AR
theta = 0.2 #MA
tau =  0.5
true_values <- c(alpha, phi, theta)
vn <- c(70, 150, 300, 500, 1000)  # Tamanhos amostrais
R <- 100 
z <- 1.96 

ar1 <- 1
ma1 <- 1

start_time <- Sys.time()

system.time({
  for (n in vn) {
    # Matrizes de resultados
    estim <- ICi <- ICs <- err <- matrix(NA, nrow = R, ncol = length(true_values))
    
    # Contadores
    calpha <- cphi <- ctheta <- 0  # Contadores para os intervalos de confiança
    bug <- 0  # Contador de falhas
    
    # Inicializa a barra de progresso
    pb <- txtProgressBar(min = 0, max = R, style = 3)
    
    for (i in 1:R) {
      y <- simu.RUR(n, phi = phi, theta = theta, alpha = alpha,tau = tau, freq = 12, link = "logit")
      fit1 <- try(suppressWarnings(RURarma.fit(y, ma = ar1, ar = ar1)), silent = TRUE)
      
      if (inherits(fit1, "try-error") ) {#|| fit1$convergence != 0) {
        bug <- bug + 1
      } else {
        estim[i, ] <- fit1$coeff
        err[i, ] <- fit1$stderror
        
        if (!any(is.na(estim[i, ])) && !any(is.na(err[i, ]))) {
          # Intervalos de Confiança
          ICi[i, ] <- estim[i, ] - (z * err[i, ])
          ICs[i, ] <- estim[i, ] + (z * err[i, ])
          
          if (ICi[i, 1] <= alpha && ICs[i, 1] >= alpha) calpha <- calpha + 1
          if (ICi[i, 2] <= phi && ICs[i, 2] >= phi) cphi <- cphi + 1
          if (ICi[i, 3] <= theta && ICs[i, 3] >= theta) ctheta <- ctheta + 1
        }
      }
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Estatísticas de desempenho
    m <- colMeans(estim, na.rm = TRUE)
    bias <- true_values - m
    biasP <- (bias / true_values) * 100
    erro <- apply(estim, 2, sd, na.rm = TRUE)
    MSE <- apply(estim, 2, var, na.rm = TRUE) + bias^2
    TC <- c(calpha, cphi, ctheta) / R
    
    # Resultados
    results <- rbind(m, bias, biasP, erro, MSE, TC)
    rownames(results) <- c("Mean", "Bias", "RB%", "SE", "MSE", "TC")
    colnames(results) <- c("alpha", "phi", "theta")
    print(c("Tamanho da Amostra:", n))
    print(round(results, 4))
    
    # Exibir avisos, se houver
    # print(warnings())
  }
})

end_time <- Sys.time()
execution_time <- end_time - start_time
print(paste("Tempo total de execução:", round(as.numeric(execution_time, units = "secs")), "segundos"))

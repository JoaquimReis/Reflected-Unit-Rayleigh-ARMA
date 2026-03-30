source("Funcoes.R")
source("SimuRUR.R")
source("FitRUR.R")

# configs gerais
set.seed(42)
n_total <- 144   # 12 anos
h       <- 12    # 1 ano
tau     <- 0.5
freq    <- 12

idx_train <- 1:(n_total - h)
idx_test  <- (n_total - h + 1):n_total

# ---------------------------------------------------------------------------
# AR(1) sem regressores
# ---------------------------------------------------------------------------

y_full  <- simu.RUR(n_total, phi = 0.8, theta = NA, alpha = 0.5,
                    tau = tau, freq = freq, link = "logit")
y_train <- ts(y_full[idx_train], frequency = freq)
y_test  <- as.numeric(y_full[idx_test])

fit1 <- RURarma.fit(y_train, ar = 1, ma = NA, tau = tau, link = "logit", h = h)

print(fit1$model)
print(round(c(loglik = fit1$loglik, AIC = fit1$aic, BIC = fit1$bic), 4))

comp1 <- data.frame(
  h        = 1:h,
  previsto = round(fit1$forecast, 4),
  real     = round(y_test, 4),
  erro     = round(fit1$forecast - y_test, 4)
)
print(comp1)


# ---------------------------------------------------------------------------
# ARMA(1,1) sem regressores
# ---------------------------------------------------------------------------

y_full2  <- simu.RUR(n_total, phi = 0.8, theta = 0.2, alpha = 0.5,
                     tau = tau, freq = freq, link = "logit")
y_train2 <- ts(y_full2[idx_train], frequency = freq)
y_test2  <- as.numeric(y_full2[idx_test])

fit2 <- RURarma.fit(y_train2, ar = 1, ma = 1, tau = tau, link = "logit", h = h)

print(fit2$model)
print(round(c(loglik = fit2$loglik, AIC = fit2$aic, BIC = fit2$bic), 4))

comp2 <- data.frame(
  h        = 1:h,
  previsto = round(fit2$forecast, 4),
  real     = round(y_test2, 4),
  erro     = round(fit2$forecast - y_test2, 4)
)
print(comp2)


# ---------------------------------------------------------------------------
# ARMA(1,1) com regressor
# ---------------------------------------------------------------------------

set.seed(99)
X_full  <- matrix(scale(cumsum(rnorm(n_total, mean = 0.02))), ncol = 1)
X_train <- X_full[idx_train, , drop = FALSE]
X_hat   <- X_full[idx_test,  , drop = FALSE]   # precisa ter exatamente h linhas

stopifnot(
  nrow(X_hat) == h,
  ncol(X_hat) == ncol(X_train)
)

fit3 <- RURarma.fit(y_train2, ar = 1, ma = 1, tau = tau, link = "logit", h = h,
                    X = X_train, X_hat = X_hat)

print(fit3$model)
print(round(c(loglik = fit3$loglik, AIC = fit3$aic, BIC = fit3$bic), 4))

comp3 <- data.frame(
  h        = 1:h,
  previsto = round(fit3$forecast, 4),
  real     = round(y_test2, 4),
  erro     = round(fit3$forecast - y_test2, 4)
)
print(comp3)


# ---------------------------------------------------------------------------
# métricas de acurácia
# ---------------------------------------------------------------------------

metricas <- function(prev, real) {
  e <- prev - real
  round(c(
    MAE  = mean(abs(e)),
    RMSE = sqrt(mean(e^2)),
    MAPE = mean(abs(e / real)) * 100
  ), 5)
}

print(rbind(
  "AR(1)"       = metricas(fit1$forecast, y_test),
  "ARMA(1,1)"   = metricas(fit2$forecast, y_test2),
  "ARMA(1,1)+X" = metricas(fit3$forecast, y_test2)
))


# ---------------------------------------------------------------------------
# diagnóstico dos resíduos quantílicos — usando fit2
# ---------------------------------------------------------------------------

resid2       <- fit2$residuals
resid2_clean <- resid2[is.finite(resid2)]

# normalidade e autocorrelação
print(shapiro.test(resid2_clean))
print(Box.test(resid2_clean, lag = 10, type = "Ljung-Box"))
print(summary(resid2_clean))


# ---------------------------------------------------------------------------
# gráficos
# ---------------------------------------------------------------------------

n_tr    <- length(y_train2)
m_order <- 1                          # max(p,q) do ARMA(1,1)
t_fit   <- (m_order + 1):n_tr
t_prev  <- (n_tr + 1):(n_tr + h)

par(mfrow = c(1, 1))

par(ask = TRUE) 


# série + ajuste + previsão
plot(1:n_tr, as.numeric(y_train2), type = "l", lwd = 1.5,
     xlim = c(1, n_tr + h), ylim = c(0, 1.3), 
     xlab = "Tempo", ylab = "y", main = "Ajuste e Previsão")
lines(t_fit, as.numeric(fit2$fitted[t_fit]), col = "blue", lwd = 1.5, lty = 2)
lines(t_prev, fit2$forecast, col = "red", lwd = 2)
points(t_prev, y_test2, pch = 16, cex = 0.8)
abline(v = n_tr, lty = 3, col = "grey50")
legend("topleft", bty = "n", lwd = c(1.5, 1.5, 2, NA), pch = c(NA, NA, NA, 16),
       col = c("black", "blue", "red", "black"), 
       cex = 0.8,
       legend = c("Observado", "Ajustado", "Previsto", "Teste"))


# resíduos ao longo do tempo
plot(resid2_clean, type = "l", lwd = 1,
     xlab = "t", ylab = "Resíduo quantílico", main = "Resíduos Quantílicos")
abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = c("red", "grey50", "red"))

# 4. qq-plot
qqnorm(resid2_clean, main = "QQ-plot", pch = 16, cex = 0.7)
qqline(resid2_clean, col = "red", lwd = 2)

par(ask = FALSE)

par(mfrow = c(1, 2))
acf(resid2_clean,  main = "ACF dos Resíduos",  lag.max = 20)
pacf(resid2_clean, main = "PACF dos Resíduos", lag.max = 20)
par(mfrow = c(1, 1))


# ---------------------------------------------------------------------------
# comparação de modelos
# ---------------------------------------------------------------------------

print(data.frame(
  modelo = c("AR(1)", "ARMA(1,1)", "ARMA(1,1)+X"),
  loglik = round(c(fit1$loglik, fit2$loglik, fit3$loglik), 4),
  k      = c(fit1$k, fit2$k, fit3$k),
  AIC    = round(c(fit1$aic, fit2$aic, fit3$aic), 4),
  BIC    = round(c(fit1$bic, fit2$bic, fit3$bic), 4),
  HQ     = round(c(fit1$hq,  fit2$hq,  fit3$hq),  4)
))


library(forecast)
library(BTSR)

source("FuncoesRUR.R")
source("FitRUR.R")

#########################
# Leitura e organizacao #
#########################

data <- read.csv("usdm_NY_none_area.csv", stringsAsFactors = FALSE)
data$date <- as.Date(data$date)
data$none <- as.numeric(data$none)
data$y_rur <- as.numeric(data$y_rur)

data <- data[is.finite(data$y_rur), ]
data$y_rur <- pmin(pmax(data$y_rur, 1e-6), 1 - 1e-6)

n <- round(nrow(data) * .8)
idx_train <- seq_len(n)
idx_test <- (n + 1):nrow(data)

Y <- data$y_rur
ytrain <- Y[idx_train]
ytest <- Y[idx_test]

date_train <- data$date[idx_train]
date_test <- data$date[idx_test]

quant <- .5
h_email <- length(ytest)
h_curto <- 12

ytrain_ts <- ts(ytrain, start = c(2000, 1), frequency = 52)

print("Tamanho das amostras")
print(c(treino = length(ytrain), teste = length(ytest)))

print("Periodo das amostras")
print(data.frame(
  amostra = c("treino", "teste"),
  inicio = c(min(date_train), min(date_test)),
  fim = c(max(date_train), max(date_test))
))

##########################################
# Covariaveis: estacoes e secas de NY    #
##########################################

# Inverno fica como categoria de referencia.
month_number <- as.integer(format(data$date, "%m"))

season <- ifelse(month_number %in% c(12, 1, 2), "inverno",
                 ifelse(month_number %in% 3:5, "primavera",
                        ifelse(month_number %in% 6:8, "verao", "outono")))

season <- factor(season,
                 levels = c("inverno", "primavera", "verao", "outono"))

season_matrix <- model.matrix(~ season)[, -1, drop = FALSE]
colnames(season_matrix) <- c("primavera", "verao", "outono")

# Janelas historicas usadas na aplicacao.
data$seca_2001_2002 <- as.numeric(
  data$date >= as.Date("2001-04-01") &
    data$date <= as.Date("2002-10-31")
)

data$seca_2016_2017 <- as.numeric(
  data$date >= as.Date("2016-06-01") &
    data$date <= as.Date("2017-02-28")
)

drought_matrix <- as.matrix(
  data[, c("seca_2001_2002", "seca_2016_2017")]
)

X0 <- cbind(season_matrix, drought_matrix)
nX <- ncol(X0)

Xtrain <- as.matrix(X0[idx_train, , drop = FALSE])
Xtest <- as.matrix(X0[idx_test, , drop = FALSE])
Xtest_curto <- Xtest[seq_len(h_curto), , drop = FALSE]

print("Covariaveis usadas no modelo")
print(colnames(X0))

print("Frequencia das dummies no treino")
print(colSums(Xtrain))

###########################
# Ordens fixas do estudo  #
###########################

# Estas ordens ficam aqui para facilitar a conferencia.
# A RUR-ARMAX(3,3) foi a que gerou os resultados citados no email.

orrurmax <- c(3, 3)
orbarmax <- c(3, 0)
orkarmax <- c(3, 0)
oruwarmax <- c(3, 2)
ormarmax <- c(3, 2)
orularmax <- c(3, 3)

# No script anterior os valores citados foram obtidos com FALSE.
# Troque para TRUE se quiser conferir o gradiente analitico.
usar_gradiente <- FALSE

##########################################
# Pequenas rotinas de apoio da aplicacao #
##########################################

lag_arg <- function(order)
{
  if(order == 0) NA else seq_len(order)
}

metricas <- function(pred, obs)
{
  pred <- as.numeric(pred)
  obs <- as.numeric(obs)

  n_eval <- min(length(pred), length(obs))

  if(n_eval == 0){
    return(data.frame(MAE = NA_real_, RMSE = NA_real_, MAPE = NA_real_))
  }

  pred <- pred[seq_len(n_eval)]
  obs <- obs[seq_len(n_eval)]

  ok <- is.finite(pred) & is.finite(obs)

  if(!any(ok)){
    return(data.frame(MAE = NA_real_, RMSE = NA_real_, MAPE = NA_real_))
  }

  obs_ok <- obs[ok]
  pred_ok <- pred[ok]

  data.frame(
    MAE = mean(abs(obs_ok - pred_ok)),
    RMSE = sqrt(mean((obs_ok - pred_ok)^2)),
    MAPE = mean(abs((obs_ok - pred_ok) / obs_ok)) * 100
  )
}

aic_btsr <- function(fit)
{
  out <- try(summary(fit)$aic, silent = TRUE)
  if(inherits(out, "try-error")) return(NA_real_)
  as.numeric(out)
}

forecast_btsr <- function(fit, h)
{
  fc <- fit$forecast

  if(is.matrix(fc) || is.data.frame(fc)){
    if("mut" %in% colnames(fc)) return(as.numeric(fc[, "mut"])[seq_len(h)])
    return(as.numeric(fc[, 1])[seq_len(h)])
  }

  as.numeric(fc)[seq_len(h)]
}

RURarma.extract <- function(yt, xreg, fit, ar, ma, tau = .5, link = "logit")
{
  yt <- as.numeric(yt)
  xreg <- as.matrix(xreg)

  stats <- make.link(link)
  linkfun <- stats$linkfun
  linkinv <- stats$linkinv

  p <- ar
  q <- ma
  n_full <- length(yt)
  m <- max(p, q)

  alpha <- as.numeric(fit$alpha)
  beta <- as.numeric(fit$beta)
  phi <- as.numeric(fit$phi)
  theta <- as.numeric(fit$theta)

  y_link <- linkfun(yt)
  Xbeta <- as.numeric(xreg %*% beta)

  eta <- rep(NA_real_, n_full)
  fitted <- rep(NA_real_, n_full)
  error <- rep(0, n_full)

  for(i in (m + 1):n_full){
    ar_part <- sum(phi * (y_link[i - seq_len(p)] - Xbeta[i - seq_len(p)]))
    ma_part <- sum(theta * error[i - seq_len(q)])

    eta[i] <- alpha + Xbeta[i] + ar_part + ma_part
    fitted[i] <- linkinv(eta[i])
    error[i] <- y_link[i] - eta[i]
  }

  list(
    mut = fitted,
    eta = eta,
    error = error,
    residuals = qnorm(pRUR(yt[(m + 1):n_full],
                           fitted[(m + 1):n_full],
                           tau = tau))
  )
}

#####################################
# Ajuste direto da RUR com regressor #
#####################################

fit_rurmax <- RURarma.fit(
  ytrain_ts,
  ar = lag_arg(orrurmax[1]),
  ma = lag_arg(orrurmax[2]),
  tau = quant,
  link = "logit",
  h = h_email,
  X = Xtrain,
  X_hat = Xtest,
  use_gradient = usar_gradiente
)

coef_rurmax <- fit_rurmax$model
beta_rows <- grep("^beta", rownames(coef_rurmax))

if(length(beta_rows) == ncol(X0)){
  rownames(coef_rurmax)[beta_rows] <- colnames(X0)
}

print("Coeficientes da RUR-ARMAX com estacoes e secas")
print(round(coef_rurmax, 5))

print("Resumo do ajuste RUR-ARMAX")
print(round(c(
  logLik = fit_rurmax$loglik,
  AIC = fit_rurmax$aic,
  BIC = fit_rurmax$bic
), 5))

##########################################
# Origem fixa curta, como teste adicional #
##########################################


fit_rurmax_h12 <- RURarma.fit(
  ytrain_ts,
  ar = lag_arg(orrurmax[1]),
  ma = lag_arg(orrurmax[2]),
  tau = quant,
  link = "logit",
  h = h_curto,
  X = Xtrain,
  X_hat = Xtest_curto,
  use_gradient = usar_gradiente
)

############################################
# Rolling/extract da RUR usando y observado #
############################################

rurmax_out <- RURarma.extract(
  yt = Y,
  xreg = X0,
  fit = fit_rurmax,
  ar = orrurmax[1],
  ma = orrurmax[2],
  tau = quant,
  link = "logit"
)

######################################
# Ajustes diretos das concorrentes   #
######################################

barmax <- BARFIMA.fit(
  yt = ytrain, xreg = Xtrain, xnew = Xtest,
  p = orbarmax[1], d = FALSE, q = orbarmax[2],
  nnew = h_email, info = TRUE, report = FALSE, linkg = "logit"
)

karmax <- KARFIMA.fit(
  yt = ytrain, xreg = Xtrain, xnew = Xtest,
  p = orkarmax[1], d = FALSE, q = orkarmax[2],
  nnew = h_email, info = TRUE, report = FALSE,
  linkg = "logit", rho = quant
)

uwarmax <- UWARFIMA.fit(
  yt = ytrain, xreg = Xtrain, xnew = Xtest,
  p = oruwarmax[1], d = FALSE, q = oruwarmax[2],
  nnew = h_email, info = TRUE, report = FALSE,
  linkg = "logit", rho = quant
)

marmax <- MARFIMA.fit(
  yt = ytrain, xreg = Xtrain, xnew = Xtest,
  p = ormarmax[1], d = FALSE, q = ormarmax[2],
  nnew = h_email, info = TRUE, report = FALSE, linkg = "logit"
)

ularmax <- ULARFIMA.fit(
  yt = ytrain, xreg = Xtrain, xnew = Xtest,
  p = orularmax[1], d = FALSE, q = orularmax[2],
  nnew = h_email, info = TRUE, report = FALSE, linkg = "logit"
)

#######################################
# Extract das concorrentes do BTSR    #
#######################################

barmax_out <- BARFIMA.extract(
  yt = Y, xreg = X0,
  coefs = list(
    alpha = barmax$coefficients[1],
    beta = barmax$coefficients[2:(nX + 1)],
    phi = barmax$coefficients[(nX + 2):(orbarmax[1] + nX + 1)],
    theta = if(orbarmax[2] == 0) {
      NULL
    } else {
      barmax$coefficients[(orbarmax[1] + nX + 2):
                            (orbarmax[1] + nX + 1 + orbarmax[2])]
    },
    nu = barmax$coefficients[(orbarmax[1] + nX + 2 + orbarmax[2])]
  )
)

karmax_out <- KARFIMA.extract(
  yt = Y, xreg = X0, rho = quant,
  coefs = list(
    alpha = karmax$coefficients[1],
    beta = karmax$coefficients[2:(nX + 1)],
    phi = karmax$coefficients[(nX + 2):(orkarmax[1] + nX + 1)],
    theta = if(orkarmax[2] == 0) {
      NULL
    } else {
      karmax$coefficients[(orkarmax[1] + nX + 2):
                            (orkarmax[1] + nX + 1 + orkarmax[2])]
    },
    nu = karmax$coefficients[(orkarmax[1] + nX + 2 + orkarmax[2])]
  )
)

uwarmax_out <- UWARFIMA.extract(
  yt = Y, xreg = X0, rho = quant,
  coefs = list(
    alpha = uwarmax$coefficients[1],
    beta = uwarmax$coefficients[2:(nX + 1)],
    phi = uwarmax$coefficients[(nX + 2):(oruwarmax[1] + nX + 1)],
    theta = if(oruwarmax[2] == 0) {
      NULL
    } else {
      uwarmax$coefficients[(oruwarmax[1] + nX + 2):
                             (oruwarmax[1] + nX + 1 + oruwarmax[2])]
    },
    nu = uwarmax$coefficients[(oruwarmax[1] + nX + 2 + oruwarmax[2])]
  )
)

#######################
# Tabelas de resultado #
#######################

resultado_aic <- data.frame(
  modelo = c("RUR-ARMAX", "BARMAX", "KARMAX", "UWARMAX",
             "MARMAX", "ULARMAX"),
  ordem = c(
    paste0("(", orrurmax[1], ",", orrurmax[2], ")"),
    paste0("(", orbarmax[1], ",", orbarmax[2], ")"),
    paste0("(", orkarmax[1], ",", orkarmax[2], ")"),
    paste0("(", oruwarmax[1], ",", oruwarmax[2], ")"),
    paste0("(", ormarmax[1], ",", ormarmax[2], ")"),
    paste0("(", orularmax[1], ",", orularmax[2], ")")
  ),
  AIC = c(
    fit_rurmax$aic,
    aic_btsr(barmax),
    aic_btsr(karmax),
    aic_btsr(uwarmax),
    aic_btsr(marmax),
    aic_btsr(ularmax)
  )
)

print("AIC dos modelos ajustados no treino")
resultado_aic_print <- resultado_aic
resultado_aic_print$AIC <- round(resultado_aic_print$AIC, 5)
print(resultado_aic_print)

resultado_origem_fixa <- rbind(
  metricas(fit_rurmax$forecast, ytest),
  metricas(forecast_btsr(barmax, h_email), ytest),
  metricas(forecast_btsr(karmax, h_email), ytest),
  metricas(forecast_btsr(uwarmax, h_email), ytest),
  metricas(forecast_btsr(marmax, h_email), ytest),
  metricas(forecast_btsr(ularmax, h_email), ytest)
)

row.names(resultado_origem_fixa) <- c(
  "RUR-ARMAX", "BARMAX", "KARMAX", "UWARMAX", "MARMAX", "ULARMAX"
)

print("Metricas fora da amostra - origem fixa para todo o teste")
print(round(resultado_origem_fixa, 4))

resultado_origem_fixa_h12 <- metricas(
  fit_rurmax_h12$forecast,
  ytest[seq_len(h_curto)]
)

row.names(resultado_origem_fixa_h12) <- "RUR-ARMAX"

print("Metricas fora da amostra - origem fixa com h = 12")
print(round(resultado_origem_fixa_h12, 4))

resultado_rolling <- rbind(
  metricas(rurmax_out$mut[idx_test], ytest),
  metricas(barmax_out$mut[idx_test], ytest),
  metricas(karmax_out$mut[idx_test], ytest),
  metricas(uwarmax_out$mut[idx_test], ytest)
)

row.names(resultado_rolling) <- c(
  "RUR-ARMAX", "BARMAX", "KARMAX", "UWARMAX"
)

print("Metricas fora da amostra - rolling/extract")
print(round(resultado_rolling, 4))


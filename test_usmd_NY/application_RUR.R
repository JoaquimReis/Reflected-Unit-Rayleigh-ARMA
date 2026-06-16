library(forecast)
library(ggplot2)
library(BTSR)

source("FuncoesRUR.R")
source("FitRUR.R")

######################
## Data preparation ##
######################

# U.S. Drought Monitor - New York
# Serie: proporcao semanal da area do estado sem seca.
# A coluna y_rur e a versao ajustada para o intervalo aberto (0,1).

data <- read.csv("4_apl_none_area.csv", stringsAsFactors = FALSE)
data$date <- as.Date(data$date)
data$none <- as.numeric(data$none)
data$y_rur <- as.numeric(data$y_rur)

data <- data[is.finite(data$y_rur), ]
data$y_rur <- pmin(pmax(data$y_rur, 1e-6), 1 - 1e-6)

n <- round(dim(data)[1] * .8)

print("Resumo da serie original: proporcao da area sem seca")
print(round(summary(data$none), 4))

print("Resumo da serie usada no modelo: y_rur em (0,1)")
print(round(summary(data$y_rur), 4))

#########################
## Train and test sets ##
#########################

datatrain <- cbind(data[1:n, ])
colnames(datatrain) <- paste0(colnames(datatrain), "_train")

datatest <- cbind(data[(n + 1):(dim(data)[1]), ])
colnames(datatest) <- paste0(colnames(datatest), "_test")

suppressMessages(attach(datatrain))
suppressMessages(attach(datatest))
suppressMessages(attach(data))

print("Tamanho das amostras")
print(c(treino = length(y_rur_train), teste = length(y_rur_test)))

######################
## Stationary tests ##
######################

truncation <- round(12 * (n / 100)^(.25)) # Schwert rule

adf.level <- tseries::adf.test(y_rur_train, k = truncation)
adf.diff <- tseries::adf.test(diff(y_rur_train), k = truncation)
kpss.level <- tseries::kpss.test(y_rur_train, lshort = FALSE)
kpss.diff <- tseries::kpss.test(diff(y_rur_train), lshort = FALSE)

table.stationarity <- data.frame(
  Series = c("In Level", "1st difference"),
  ADF = c(adf.level$statistic, adf.diff$statistic),
  p.value.ADF = c(adf.level$p.value, adf.diff$p.value),
  KPSS = c(kpss.level$statistic, kpss.diff$statistic),
  p.value.KPSS = c(kpss.level$p.value, kpss.diff$p.value)
)

print("Testes de estacionariedade")
print(data.frame(
  Series = table.stationarity$Series,
  round(table.stationarity[, -1], 4)
))

# ACF and PACF of the time series
acf(y_rur_train, main = "ACF - serie em nivel")
pacf(y_rur_train, main = "PACF - serie em nivel")
acf(diff(y_rur_train), main = "ACF - primeira diferenca")
pacf(diff(y_rur_train), main = "PACF - primeira diferenca")

########################
## Auxiliary routines ##
########################

safe_aic <- function(fit)
{
  if(inherits(fit, "try-error") || is.null(fit)) return(Inf)
  sm <- try(summary(fit), silent = TRUE)
  if(inherits(sm, "try-error") || is.null(sm$aic) || !is.finite(sm$aic)) return(Inf)
  sm$aic
}

fit_ok <- function(fit)
{
  !inherits(fit, "try-error") &&
    !is.null(fit$convergence) &&
    fit$convergence == 0 &&
    is.finite(safe_aic(fit))
}

forecast_mean <- function(fit, h)
{
  fc <- fit$forecast
  if(is.null(fc)) return(rep(NA_real_, h))

  if(is.matrix(fc) || is.data.frame(fc)){
    if("mut" %in% colnames(fc)) return(as.numeric(fc[, "mut"])[1:h])
    return(as.numeric(fc[, 1])[1:h])
  }

  as.numeric(fc)[1:h]
}

coef_btsr <- function(fit)
{
  ans <- try(t(summary(fit)$coefficients[, c(1, 4), drop = FALSE]), silent = TRUE)
  if(inherits(ans, "try-error")) return(matrix(NA, nrow = 0, ncol = 0))
  ans
}

coef_rur <- function(fit)
{
  ans <- try(t(fit$model[, c(1, 4), drop = FALSE]), silent = TRUE)
  if(inherits(ans, "try-error")) return(matrix(NA, nrow = 0, ncol = 0))
  ans
}

coef_arima <- function(fit)
{
  ans <- try(lmtest::coeftest(fit), silent = TRUE)
  if(inherits(ans, "try-error") || length(ans) == 0){
    return(matrix(NA, nrow = 0, ncol = 0))
  }

  try_ans <- try(t(ans[, c(1, 4), drop = FALSE]), silent = TRUE)
  if(inherits(try_ans, "try-error")) return(matrix(NA, nrow = 0, ncol = 0))
  try_ans
}

lag_arg <- function(order)
{
  if(order == 0) NA else 1:order
}

########################
## Fitting the models ##
########################

quant <- .5
h <- length(y_rur_test)

# Fitting ARIMA and SARIMA as Gaussian benchmarks
a01 <- auto.arima(y_rur_train, allowdrift = FALSE)
new1 <- forecast(a01, h = h)

s01 <- auto.arima(ts(y_rur_train, frequency = 52), allowdrift = FALSE)
snew1 <- forecast(s01, h = h)

# Fitting the unit ARMA models
order <- matrix(NA, 15, 7)
colnames(order) <- c("p", "q", "BARMA", "KARMA", "UWARMA",
                     "MARMA", "ULARMA")

cont <- 1

for(i in 0:3){
  for(j in 0:3){
    if(i == 0 && j == 0) next

    barma <- try(BARFIMA.fit(yt = y_rur_train, p = i, d = FALSE, q = j,
                             nnew = h, info = TRUE, report = FALSE),
                 silent = TRUE)

    karma <- try(suppressWarnings(KARFIMA.fit(yt = y_rur_train, p = i, d = FALSE, q = j,
                                              nnew = h, info = TRUE, rho = quant,
                                              control = list(method = "Nelder-Mead",
                                                             stopcr = 1e-2),
                                              report = FALSE)),
                 silent = TRUE)

    uwarma <- try(suppressWarnings(UWARFIMA.fit(yt = y_rur_train, p = i, d = FALSE, q = j,
                                                nnew = h, info = TRUE, rho = quant,
                                                report = FALSE)),
                  silent = TRUE)

    marma <- try(suppressWarnings(MARFIMA.fit(yt = y_rur_train, p = i, d = FALSE, q = j,
                                              nnew = h, info = TRUE, report = FALSE)),
                 silent = TRUE)

    ularma <- try(suppressWarnings(ULARFIMA.fit(yt = y_rur_train, p = i, d = FALSE, q = j,
                                                nnew = h, info = TRUE, report = FALSE)),
                  silent = TRUE)

    order[cont, ] <- c(i, j,
                       if(fit_ok(barma)) safe_aic(barma) else Inf,
                       if(fit_ok(karma)) safe_aic(karma) else Inf,
                       if(fit_ok(uwarma)) safe_aic(uwarma) else Inf,
                       if(fit_ok(marma)) safe_aic(marma) else Inf,
                       if(fit_ok(ularma)) safe_aic(ularma) else Inf)
    cont <- cont + 1
  }
}

# RUR is fitted separately because it is our implementation
order_rur <- matrix(NA, 15, 3)
colnames(order_rur) <- c("p", "q", "RUR")

cont <- 1
for(i in 0:3){
  for(j in 0:3){
    if(i == 0 && j == 0) next

    rur_tmp <- try(suppressWarnings(RURarma.fit(
      ts(y_rur_train, start = c(2000, 1), frequency = 52),
      ar = lag_arg(i),
      ma = lag_arg(j),
      tau = quant,
      link = "logit",
      h = h,
      use_gradient = FALSE
    )), silent = TRUE)

    order_rur[cont, ] <- c(i, j,
                           if(!inherits(rur_tmp, "try-error") &&
                              !is.null(rur_tmp$conv) &&
                              rur_tmp$conv == 0 &&
                              is.finite(rur_tmp$aic)) rur_tmp$aic else Inf)
    cont <- cont + 1
  }
}

print("Grade de AIC para modelos unitarios do BTSR")
print(round(order, 4))

print("Grade de AIC para modelos RUR")
print(round(order_rur, 4))

#######################################
## Selecting the order of each model ##
#######################################

orbarma <- order[which(order[, "BARMA"] == min(order[, "BARMA"], na.rm = TRUE))[1], c("p", "q")]
orkarma <- order[which(order[, "KARMA"] == min(order[, "KARMA"], na.rm = TRUE))[1], c("p", "q")]
oruwarma <- order[which(order[, "UWARMA"] == min(order[, "UWARMA"], na.rm = TRUE))[1], c("p", "q")]
ormarma <- order[which(order[, "MARMA"] == min(order[, "MARMA"], na.rm = TRUE))[1], c("p", "q")]
orularma <- order[which(order[, "ULARMA"] == min(order[, "ULARMA"], na.rm = TRUE))[1], c("p", "q")]
orrur <- order_rur[which(order_rur[, "RUR"] == min(order_rur[, "RUR"], na.rm = TRUE))[1], c("p", "q")]

print("Ordens selecionadas por AIC")
print(rbind(
  BARMA = orbarma,
  KARMA = orkarma,
  UWARMA = oruwarma,
  MARMA = ormarma,
  ULARMA = orularma,
  RUR = orrur
))

##########################
## Final model fittings ##
##########################

barma <- BARFIMA.fit(yt = y_rur_train, p = orbarma[1], d = FALSE, q = orbarma[2],
                     nnew = h, info = TRUE, report = FALSE)

karma <- KARFIMA.fit(yt = y_rur_train, p = orkarma[1], d = FALSE, q = orkarma[2],
                     nnew = h, rho = quant,
                     control = list(method = "Nelder-Mead", stopcr = 1e-2),
                     info = TRUE, report = FALSE)

uwarma <- UWARFIMA.fit(yt = y_rur_train, p = oruwarma[1], d = FALSE, q = oruwarma[2],
                       nnew = h, rho = quant, info = TRUE, report = FALSE)

marma <- MARFIMA.fit(yt = y_rur_train, p = ormarma[1], d = FALSE, q = ormarma[2],
                     nnew = h, info = TRUE, report = FALSE)

ularma <- ULARFIMA.fit(yt = y_rur_train, p = orularma[1], d = FALSE, q = orularma[2],
                       nnew = h, info = TRUE, report = FALSE)

rurarma <- RURarma.fit(ts(y_rur_train, start = c(2000, 1), frequency = 52),
                       ar = lag_arg(orrur[1]),
                       ma = lag_arg(orrur[2]),
                       tau = quant,
                       link = "logit",
                       h = h,
                       use_gradient = FALSE)

##################
## Coefficients ##
##################

barma_coeff <- coef_btsr(barma)
karma_coeff <- coef_btsr(karma)
uwarma_coeff <- coef_btsr(uwarma)
marma_coeff <- coef_btsr(marma)
ularma_coeff <- coef_btsr(ularma)
rur_coeff <- coef_rur(rurarma)
sarima_coeff <- coef_arima(s01)
arima_coeff <- coef_arima(a01)

print("Coeficientes BARMA")
print(round(barma_coeff, 4))

print("Coeficientes KARMA")
print(round(karma_coeff, 4))

print("Coeficientes UWARMA")
print(round(uwarma_coeff, 4))

print("Coeficientes MARMA")
print(round(marma_coeff, 4))

print("Coeficientes ULARMA")
print(round(ularma_coeff, 4))

print("Coeficientes RUR-ARMA")
print(round(rur_coeff, 4))

print("Coeficientes SARIMA")
print(round(sarima_coeff, 4))

print("Coeficientes ARIMA")
print(round(arima_coeff, 4))

##########################
## Residual diagnostics ##
##########################

# Ljung-Box test. p-values altos indicam menor evidencia de autocorrelacao residual.
print("Ljung-Box residual tests")
print(Box.test(barma$residuals, lag = 20, fitdf = sum(orbarma)))
print(Box.test(karma$residuals, lag = 20, fitdf = sum(orkarma)))
print(Box.test(uwarma$residuals, lag = 20, fitdf = sum(oruwarma)))
print(Box.test(marma$residuals, lag = 20, fitdf = sum(ormarma)))
print(Box.test(ularma$residuals, lag = 20, fitdf = sum(orularma)))
print(Box.test(rurarma$residuals, lag = 20, fitdf = sum(orrur)))
print(Box.test(s01$residuals, lag = 20, fitdf = length(coef(s01))))

# ACF and PACF of the residuals
acf(rurarma$residuals, main = "ACF - residuos RUR")
pacf(rurarma$residuals, main = "PACF - residuos RUR")
acf(uwarma$residuals, main = "ACF - residuos Unit-Weibull")
pacf(uwarma$residuals, main = "PACF - residuos Unit-Weibull")
acf(barma$residuals, main = "ACF - residuos Beta")
pacf(barma$residuals, main = "PACF - residuos Beta")

#########################
## Out-of-sample error ##
#########################

pred_barma <- forecast_mean(barma, h)
pred_karma <- forecast_mean(karma, h)
pred_uwarma <- forecast_mean(uwarma, h)
pred_marma <- forecast_mean(marma, h)
pred_ularma <- forecast_mean(ularma, h)
pred_rur <- as.numeric(rurarma$forecast)

results_outsample <- rbind(
  forecast::accuracy(pred_barma, y_rur_test),
  forecast::accuracy(pred_karma, y_rur_test),
  forecast::accuracy(pred_uwarma, y_rur_test),
  forecast::accuracy(pred_marma, y_rur_test),
  forecast::accuracy(pred_ularma, y_rur_test),
  forecast::accuracy(pred_rur, y_rur_test),
  forecast::accuracy(as.numeric(new1$mean), y_rur_test),
  forecast::accuracy(as.numeric(snew1$mean), y_rur_test)
)[, c(3, 2, 5)]

row.names(results_outsample) <- c("BARMA", "KARMA", "UWARMA",
                                  "MARMA", "ULARMA", "RUR-ARMA",
                                  "ARIMA", "SARIMA")

print("Erros fora da amostra")
print(round(results_outsample, 4))

###########################
## AIC comparison table  ##
###########################

table_aic <- data.frame(
  model = c("BARMA", "KARMA", "UWARMA", "MARMA", "ULARMA",
            "RUR-ARMA", "ARIMA", "SARIMA"),
  AIC = c(
    safe_aic(barma),
    safe_aic(karma),
    safe_aic(uwarma),
    safe_aic(marma),
    safe_aic(ularma),
    rurarma$aic,
    AIC(a01),
    AIC(s01)
  ),
  MAE = results_outsample[, "MAE"],
  RMSE = results_outsample[, "RMSE"],
  MAPE = results_outsample[, "MAPE"]
)

table_aic <- table_aic[order(table_aic$AIC), ]

print("Tabela final ordenada por AIC")
table_aic_print <- table_aic
table_aic_print[, c("AIC", "MAE", "RMSE", "MAPE")] <-
  round(table_aic_print[, c("AIC", "MAE", "RMSE", "MAPE")], 4)
print(table_aic_print)

# Nesta aplicacao, a RUR tem menor AIC na grade considerada.
# O MAPE fica instavel porque existem valores muito proximos de zero.

############################
## Comparison bar chart   ##
############################

df_aic <- table_aic
df_aic$model <- factor(df_aic$model, levels = rev(df_aic$model))

ggplot(df_aic, aes(x = model, y = AIC, fill = model)) +
  geom_col(width = .75) +
  coord_flip() +
  labs(y = "AIC", x = "", fill = "") +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey85"),
    axis.title.y = element_text(face = "bold", color = "black", size = 9),
    axis.title.x = element_text(face = "bold", color = "black", size = 9),
    axis.text.x = element_text(face = "bold", color = "black", size = 8),
    axis.text.y = element_text(face = "bold", color = "black", size = 8)
  )

###########################
## RUR fitted values plot ##
###########################

fit_rur_data <- data.frame(
  date = date_train,
  Observado = y_rur_train,
  RUR = as.numeric(rurarma$fitted)
)

fit_rur_data <- fit_rur_data[is.finite(fit_rur_data$RUR), ]

fit_rur_long <- data.frame(
  date = rep(fit_rur_data$date, 2),
  value = c(fit_rur_data$Observado, fit_rur_data$RUR),
  Serie = rep(c("Observado", "RUR ajustada"), each = nrow(fit_rur_data))
)

ggplot(fit_rur_long, aes(x = date, y = value, color = Serie)) +
  geom_line(linewidth = .55) +
  labs(title = "Ajuste da RUR",
       y = "Proporcao", x = "", color = "") +
  scale_color_manual(values = c("Observado" = "grey25",
                                "RUR ajustada" = "#0072B2")) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(colour = "grey90"),
    plot.title = element_text(face = "bold", color = "black", size = 12),
    plot.subtitle = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold", color = "black", size = 9),
    axis.title.x = element_text(face = "bold", color = "black", size = 9),
    axis.text.x = element_text(face = "bold", color = "black", size = 8),
    axis.text.y = element_text(face = "bold", color = "black", size = 8),
    legend.text = element_text(face = "bold", size = 8)
  )

###########################
## RUR fitted plot: zoom ##
###########################

last_weeks <- 5 * 52
fit_rur_zoom <- tail(fit_rur_data, last_weeks)

fit_rur_zoom_long <- data.frame(
  date = rep(fit_rur_zoom$date, 2),
  value = c(fit_rur_zoom$Observado, fit_rur_zoom$RUR),
  Serie = rep(c("Observado", "RUR ajustada"), each = nrow(fit_rur_zoom))
)

ggplot(fit_rur_zoom_long, aes(x = date, y = value, color = Serie)) +
  geom_line(linewidth = .7) +
  labs(title = "Ajuste da RUR no final do treino",
       y = "Proporcao", x = "", color = "") +
  scale_color_manual(values = c("Observado" = "grey25",
                                "RUR ajustada" = "#0072B2")) +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(colour = "grey90"),
    plot.title = element_text(face = "bold", color = "black", size = 12),
    plot.subtitle = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold", color = "black", size = 9),
    axis.title.x = element_text(face = "bold", color = "black", size = 9),
    axis.text.x = element_text(face = "bold", color = "black", size = 8),
    axis.text.y = element_text(face = "bold", color = "black", size = 8),
    legend.text = element_text(face = "bold", size = 8)
  )

############################
## Test-set visualization ##
############################

plot_data <- data.frame(
  date = date_test,
  Observed = y_rur_test,
  RUR = pred_rur,
  UWARMA = pred_uwarma,
  Matsuoka = pred_marma
)

plot_data_long <- data.frame(
  date = rep(plot_data$date, 4),
  value = c(plot_data$Observed, plot_data$RUR,
            plot_data$UWARMA, plot_data$Matsuoka),
  model = rep(c("Observed", "RUR", "UWARMA", "Matsuoka"),
              each = nrow(plot_data))
)

ggplot(plot_data_long, aes(x = date, y = value, color = model)) +
  geom_line(linewidth = .6) +
  labs(title = "Previsoes no conjunto de teste",
       y = "Proporcao", x = "", color = "") +
  theme(
    legend.position = "bottom",
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_line(colour = "grey90"),
    plot.title = element_text(face = "bold", color = "black", size = 12),
    plot.subtitle = element_text(color = "black", size = 10),
    axis.title.y = element_text(face = "bold", color = "black", size = 9),
    axis.title.x = element_text(face = "bold", color = "black", size = 9),
    axis.text.x = element_text(face = "bold", color = "black", size = 8),
    axis.text.y = element_text(face = "bold", color = "black", size = 8),
    legend.text = element_text(face = "bold", size = 8)
  )

# The fixed-origin forecasts tend to flatten over the test period.
# This plot is useful for prediction behavior, not for showing one-step fitted values.

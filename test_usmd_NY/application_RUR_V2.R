library(forecast)
library(ggplot2)

source("FuncoesRUR.R")
source("FitRUR.R")

######################
## Data preparation ##
######################

# U.S. Drought Monitor - New York
# Serie semanal: proporcao da area do estado sem seca.

data <- read.csv("usdm_NY_none_area.csv", stringsAsFactors = FALSE)
data$date <- as.Date(data$date)
data$none <- as.numeric(data$none)
data$y_rur <- as.numeric(data$y_rur)

data <- data[is.finite(data$y_rur), ]
data$y_rur <- pmin(pmax(data$y_rur, 1e-6), 1 - 1e-6)

n <- round(nrow(data) * .8)
idx_train <- seq_len(n)
idx_test <- (n + 1):nrow(data)

ytrain <- data$y_rur[idx_train]
ytest <- data$y_rur[idx_test]
date_train <- data$date[idx_train]
date_test <- data$date[idx_test]
h <- length(ytest)

ytrain_ts <- ts(ytrain, start = c(2000, 1), frequency = 52)

quant <- .5
pmax <- 3
qmax <- 3

print("Tamanho das amostras")
print(c(treino = length(ytrain), teste = length(ytest)))

#################################
## Seasonal dummy construction ##
#################################

# Inverno e a categoria de referencia.
# Os coeficientes de primavera, verao e outono devem ser interpretados
# em comparacao com o inverno, mantendo a estrutura ARMA constante.

month_number <- as.integer(format(data$date, "%m"))

season <- ifelse(month_number %in% c(12, 1, 2), "inverno",
                 ifelse(month_number %in% 3:5, "primavera",
                        ifelse(month_number %in% 6:8, "verao", "outono")))

season <- factor(season,
                 levels = c("inverno", "primavera", "verao", "outono"))

season_matrix <- model.matrix(~ season)[, -1, drop = FALSE]
colnames(season_matrix) <- c("primavera", "verao", "outono")

################################
## Historical drought dummies ##
################################

# As duas janelas representam as grandes quedas discutidas no estudo 2.0.
# Elas sao intervencoes historicas e nao previsores de novas secas.

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

##################################
## Covariate scenarios          ##
##################################

# O quarto cenario e a uniao das duas melhores alternativas do estudo 2.0.

covariate_sets <- list(
  sem_covariaveis = NULL,
  estacoes = season_matrix,
  secas_historicas = drought_matrix,
  estacoes_e_secas = cbind(season_matrix, drought_matrix)
)

########################
## Auxiliary routines ##
########################

lag_arg <- function(order)
{
  if(order == 0) NA else seq_len(order)
}

accuracy_row <- function(pred, obs)
{
  pred <- as.numeric(pred)
  obs <- as.numeric(obs)
  ok <- is.finite(pred) & is.finite(obs)

  if(!any(ok)) return(c(MAE = Inf, RMSE = Inf))

  c(
    MAE = mean(abs(obs[ok] - pred[ok])),
    RMSE = sqrt(mean((obs[ok] - pred[ok])^2))
  )
}

residual_diagnostics <- function(residuals, fitdf = 0, lag_max = 52)
{
  residuals <- as.numeric(residuals)
  residuals <- residuals[is.finite(residuals)]

  lb_value <- function(lag)
  {
    adjusted_df <- min(fitdf, lag - 1)
    test <- Box.test(residuals, lag = lag, type = "Ljung-Box",
                     fitdf = adjusted_df)
    as.numeric(test$p.value)
  }

  acf_values <- as.numeric(
    acf(residuals, lag.max = lag_max, plot = FALSE)$acf[-1]
  )

  band <- 1.96 / sqrt(length(residuals))

  c(
    LB20 = lb_value(20),
    LB52 = lb_value(52),
    max_ACF = max(abs(acf_values)),
    lags_fora = sum(abs(acf_values) > band)
  )
}

#################################
## RUR fitting for one scenario ##
#################################

fit_rur_once <- function(p, q, scenario, Xall = NULL)
{
  if(is.null(Xall)){
    Xtrain <- NULL
    Xtest <- NULL
    k <- 0
  } else {
    Xtrain <- as.matrix(Xall[idx_train, , drop = FALSE])
    Xtest <- as.matrix(Xall[idx_test, , drop = FALSE])
    k <- ncol(Xtrain)
  }

  fit <- try(suppressWarnings(RURarma.fit(
    ytrain_ts,
    ar = lag_arg(p),
    ma = lag_arg(q),
    tau = quant,
    link = "logit",
    h = h,
    X = if(is.null(Xtrain)) NA else Xtrain,
    X_hat = if(is.null(Xtest)) NA else Xtest,
    use_gradient = FALSE
  )), silent = TRUE)

  ok <- !inherits(fit, "try-error") &&
    !is.null(fit$conv) && fit$conv == 0 &&
    !is.null(fit$aic) && is.finite(fit$aic)

  if(!ok){
    row <- data.frame(
      cenario = scenario, p = p, q = q, k = k,
      logLik = NA, AIC = Inf, BIC = Inf,
      MAE = Inf, RMSE = Inf,
      LB20 = NA, LB52 = NA, max_ACF = NA, lags_fora = NA,
      convergiu = FALSE
    )
    return(list(fit = fit, row = row))
  }

  errors <- accuracy_row(fit$forecast, ytest)
  diag <- residual_diagnostics(fit$residuals, fitdf = p + q + k)

  row <- data.frame(
    cenario = scenario, p = p, q = q, k = k,
    logLik = fit$loglik, AIC = fit$aic, BIC = fit$bic,
    MAE = unname(errors["MAE"]),
    RMSE = unname(errors["RMSE"]),
    LB20 = unname(diag["LB20"]),
    LB52 = unname(diag["LB52"]),
    max_ACF = unname(diag["max_ACF"]),
    lags_fora = unname(diag["lags_fora"]),
    convergiu = TRUE
  )

  list(fit = fit, row = row)
}

fit_rur_grid <- function(scenario, Xall = NULL)
{
  rows <- list()
  fits <- list()
  cont <- 1

  for(p in 0:pmax){
    for(q in 0:qmax){
      if(p == 0 && q == 0) next

      ans <- fit_rur_once(p, q, scenario, Xall)
      rows[[cont]] <- ans$row
      fits[[cont]] <- ans$fit
      cont <- cont + 1
    }
  }

  grid <- do.call(rbind, rows)
  valid <- which(is.finite(grid$AIC))

  if(length(valid) == 0){
    return(list(fit = NULL, summary = grid[1, ], grid = grid))
  }

  best <- valid[which.min(grid$AIC[valid])]

  list(
    fit = fits[[best]],
    summary = grid[best, ],
    grid = grid[order(grid$AIC), ]
  )
}

####################################
## Complete order search          ##
####################################

# Cada cenario recebe sua propria busca de p e q.

results <- list()

for(scenario in names(covariate_sets)){
  print(paste("Ajustando a grade completa para:", scenario))
  results[[scenario]] <- fit_rur_grid(
    scenario, covariate_sets[[scenario]]
  )
}

comparison <- do.call(rbind, lapply(results, function(x) x$summary))
comparison <- comparison[order(comparison$AIC), ]

comparison_print <- comparison
numeric_columns <- c("logLik", "AIC", "BIC", "MAE", "RMSE",
                     "LB20", "LB52", "max_ACF")
comparison_print[, numeric_columns] <-
  round(comparison_print[, numeric_columns], 5)

print("Comparacao dos quatro cenarios")
print(comparison_print)

####################################
## Combined model coefficients    ##
####################################

combined_result <- results[["estacoes_e_secas"]]
combined_fit <- combined_result$fit
combined_X <- covariate_sets[["estacoes_e_secas"]]

combined_coefficients <- combined_fit$model
beta_rows <- grep("^beta", rownames(combined_coefficients))

if(length(beta_rows) == ncol(combined_X)){
  rownames(combined_coefficients)[beta_rows] <- colnames(combined_X)
}

print("Coeficientes da RUR com estacoes e secas")
print(round(combined_coefficients, 5))

# Para as dummies de seca, um coeficiente negativo indica menor proporcao
# condicional sem seca. Um sinal positivo pede cautela: a janela pode incluir
# semanas de recuperacao ou seu efeito pode estar sendo absorvido pelos lags.

########################################
## Likelihood-ratio tests              ##
########################################

# Os modelos abaixo precisam ter a mesma ordem ARMA. Como a ordem (3,3) do
# modelo combinado pode nao convergir nos modelos reduzidos, colquei a
# melhor ordem que tenha convergido nos quatro cenarios.

valid_order_keys <- lapply(results, function(x){
  valid <- x$grid[is.finite(x$grid$AIC), c("p", "q")]
  paste(valid$p, valid$q, sep = "_")
})

common_order_keys <- Reduce(intersect, valid_order_keys)

if(length(common_order_keys) == 0){
  stop("Nao foi encontrada uma ordem comum com convergencia nos quatro cenarios.")
}

combined_grid <- combined_result$grid
combined_grid$order_key <- paste(combined_grid$p, combined_grid$q, sep = "_")
combined_common <- combined_grid[
  combined_grid$order_key %in% common_order_keys & is.finite(combined_grid$AIC),
]
combined_common <- combined_common[order(combined_common$AIC), ]

common_p <- combined_common$p[1]
common_q <- combined_common$q[1]

print("Ordem comum usada nos testes de razao de verossimilhancas")
print(c(p = common_p, q = common_q))

same_order_results <- list()

for(scenario in names(covariate_sets)){
  same_order_results[[scenario]] <- fit_rur_once(
    common_p, common_q, scenario, covariate_sets[[scenario]]
  )
}

lr_test <- function(reduced_name, complete_name, degrees_freedom)
{
  reduced <- same_order_results[[reduced_name]]$fit
  complete <- same_order_results[[complete_name]]$fit

  if(inherits(reduced, "try-error") || inherits(complete, "try-error") ||
     is.null(reduced$loglik) || is.null(complete$loglik)){
    return(data.frame(
      comparacao = paste(reduced_name, "vs", complete_name),
      LR = NA, graus_liberdade = degrees_freedom, p_value = NA
    ))
  }

  statistic <- 2 * (complete$loglik - reduced$loglik)

  data.frame(
    comparacao = paste(reduced_name, "vs", complete_name),
    LR = statistic,
    graus_liberdade = degrees_freedom,
    p_value = pchisq(statistic, df = degrees_freedom, lower.tail = FALSE)
  )
}

lr_results <- rbind(
  lr_test("sem_covariaveis", "estacoes_e_secas", 5),
  lr_test("estacoes", "estacoes_e_secas", 2),
  lr_test("secas_historicas", "estacoes_e_secas", 3)
)

lr_results$LR <- round(lr_results$LR, 5)
lr_results$p_value <- round(lr_results$p_value, 6)

print("Testes de razao de verossimilhancas na ordem do modelo combinado")
print(lr_results)

# No teste estacoes vs estacoes_e_secas, p < 0.05 as duas dummies
# de seca acrescentam informacao

####################################
## Residual diagnostics           ##
####################################

residual_table <- comparison[, c("cenario", "p", "q", "LB20", "LB52",
                                  "max_ACF", "lags_fora")]
residual_table <- residual_table[order(residual_table$lags_fora), ]
residual_table[, c("LB20", "LB52", "max_ACF")] <-
  round(residual_table[, c("LB20", "LB52", "max_ACF")], 6)

print("Comparacao dos residuos quantilicos")
print(residual_table)

par(mfrow = c(2, 2))
acf(combined_fit$residuals, lag.max = 52,
    main = "ACF - RUR com estacoes e secas")
pacf(combined_fit$residuals, lag.max = 52,
     main = "PACF - RUR com estacoes e secas")
qqnorm(combined_fit$residuals,
       main = "Normal Q-Q - residuos quantilicos")
qqline(combined_fit$residuals, col = "red")
hist(combined_fit$residuals, breaks = 25,
     main = "Histograma - residuos quantilicos", xlab = "Residuo")
par(mfrow = c(1, 1))

####################################
## AIC graphic     ##
####################################

g_aic <- ggplot(comparison,
                aes(x = reorder(cenario, AIC), y = AIC, fill = cenario)) +
  geom_col(width = .72) +
  coord_flip() +
  labs(title = "RUR: estacoes e dummies de seca",
       subtitle = "Comparacao dos quatro cenarios por AIC",
       x = "", y = "AIC") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

print(g_aic)

####################################
## Fitted values comparison       ##
####################################

baseline_fit <- results[["sem_covariaveis"]]$fit
season_fit <- results[["estacoes"]]$fit
drought_fit <- results[["secas_historicas"]]$fit

fitted_data <- rbind(
  data.frame(date = date_train, valor = ytrain, serie = "Observado"),
  data.frame(date = date_train, valor = as.numeric(baseline_fit$fitted),
             serie = "RUR original"),
  data.frame(date = date_train, valor = as.numeric(combined_fit$fitted),
             serie = "RUR com estacoes e secas")
)

fitted_data <- fitted_data[is.finite(fitted_data$valor), ]

g_fitted <- ggplot(fitted_data, aes(x = date, y = valor, color = serie)) +
  geom_line(linewidth = .52) +
  scale_color_manual(values = c(
    "Observado" = "grey25",
    "RUR original" = "#0072B2",
    "RUR com estacoes e secas" = "#D55E00"
  )) +
  labs(title = "Ajuste da RUR dentro da amostra",
       x = "", y = "Proporcao", color = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(g_fitted)

####################################
## Zoom around historical droughts ##
####################################

zoom_data <- fitted_data[
  fitted_data$date >= as.Date("2015-01-01") &
    fitted_data$date <= as.Date("2018-01-01"),
]

g_zoom <- ggplot(zoom_data, aes(x = date, y = valor, color = serie)) +
  annotate("rect", xmin = as.Date("2016-06-01"),
           xmax = as.Date("2017-02-28"),
           ymin = -Inf, ymax = Inf, fill = "#D55E00", alpha = .10) +
  geom_line(linewidth = .7) +
  scale_color_manual(values = c(
    "Observado" = "grey25",
    "RUR original" = "#0072B2",
    "RUR com estacoes e secas" = "#D55E00"
  )) +
  labs(title = "Ajuste durante a seca de 2016-2017",
       subtitle = "A faixa sombreada representa a janela da dummy de seca",
       x = "", y = "Proporcao", color = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(g_zoom)

####################################
## Test-set forecasts             ##
####################################

prediction_data <- rbind(
  data.frame(date = date_test, valor = ytest, serie = "Observado"),
  data.frame(date = date_test, valor = as.numeric(baseline_fit$forecast),
             serie = "RUR original"),
  data.frame(date = date_test, valor = as.numeric(season_fit$forecast),
             serie = "RUR com estacoes"),
  data.frame(date = date_test, valor = as.numeric(drought_fit$forecast),
             serie = "RUR com secas"),
  data.frame(date = date_test, valor = as.numeric(combined_fit$forecast),
             serie = "RUR com estacoes e secas")
)

g_prediction <- ggplot(prediction_data,
                       aes(x = date, y = valor, color = serie)) +
  geom_line(linewidth = .58) +
  labs(title = "Previsoes no conjunto de teste",
       subtitle = "Previsoes de origem fixa para todo o horizonte",
       x = "", y = "Proporcao", color = "") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold"))

print(g_prediction)





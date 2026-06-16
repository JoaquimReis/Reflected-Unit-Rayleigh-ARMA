# ---------------------------------------------------------------------------
# Funcoes da distribuicao Reflected Unit Rayleigh (RUR)
# ---------------------------------------------------------------------------

validate_RUR_args <- function(x = NULL, mu = 0.7, tau = 0.5)
{
  if(any(mu <= 0 | mu >= 1, na.rm = TRUE)){
    stop("mu must be between 0 and 1")
  }

  if(any(tau <= 0 | tau >= 1, na.rm = TRUE)){
    stop("tau must be between 0 and 1")
  }

  if(!is.null(x) && any(x <= 0 | x >= 1, na.rm = TRUE)){
    stop("x must be between 0 and 1")
  }
}

# Density function
dRUR <- function(y, mu = 0.7, tau = 0.5, log = FALSE)
{
  validate_RUR_args(x = y, mu = mu, tau = tau)

  log_fy <- log(2) + log(-log(1 - tau)) - log(1 - y) -
    2 * log(-log(1 - mu)) + log(-log(1 - y)) +
    log(1 - tau) * ((log(1 - y) / log(1 - mu))^2)

  if(log) log_fy else exp(log_fy)
}

# Cumulative distribution function
pRUR <- function(q, mu = 0.7, tau = 0.5,
                 lower.tail = TRUE, log.p = FALSE)
{
  validate_RUR_args(x = q, mu = mu, tau = tau)

  cdf <- 1 - (1 - tau)^((log(1 - q) / log(1 - mu))^2)
  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) log(cdf) else cdf
}

# Quantile function
qRUR <- function(u, mu = 0.7, tau = 0.5,
                 lower.tail = TRUE, log.p = FALSE)
{
  if(log.p) u <- exp(u)
  if(!lower.tail) u <- 1 - u

  validate_RUR_args(x = u, mu = mu, tau = tau)

  exponent <- sqrt(log(1 - u) / log(1 - tau))
  1 - (1 - mu)^exponent
}

# Random generation by inverse transform
rRUR <- function(n, mu = 0.7, tau = 0.5)
{
  if(length(n) != 1 || n < 0) stop("n must be a non-negative scalar")
  validate_RUR_args(mu = mu, tau = tau)

  y <- qRUR(runif(n), mu = mu, tau = tau)
  pmin(pmax(y, 1e-6), 1 - 1e-6)
}

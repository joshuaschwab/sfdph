rexpit <- function(x) {
  y <- rep(NA_integer_, length(x))
  index <- !is.na(x)
  y[index] <- rbinom(n = sum(index), size = 1, prob = plogis(x[index]))
  return(y)
}

Bound <- function(x, bounds) {
  x[x < min(bounds)] <- min(bounds)
  x[x > max(bounds)] <- max(bounds)
  return(x)
}

GenerateData <- function(n, num.time.points, censoring, repeated.y) {
  sex <- rbinom(n, size = 1, prob = 0.5)
  age <- round(runif(n, min = 20, max = 80))
  prev.renal <- rbinom(n, size = 1, prob = 0.3)
  prev.ldl <- Bound(round(rnorm(n, 120, 25)), c(50, 500))
  prev.abc <- y <- 0
  censored <- rep(F, n)
  df <- data.frame(sex, age, renal0 = prev.renal, ldl0 = prev.ldl)
  index <- rep(T, n) #index of uncensored and alive patients
  for (t in 1:num.time.points) {
    renal <- ldl <- abc <- rep(NA, n)
    renal[index] <- renal[index] <- rexpit(-2 + 0.3 * sex + 0.01 * age + prev.renal)[index]
    ldl[index] <- round((prev.ldl + rnorm(n, mean = 0, sd = 10) + 15 * prev.abc)[index])
    abc[index] <- rexpit(1.1 * renal + 0.0005 * (ldl - 120)^2)[index]
    temp.df <- data.frame(renal, ldl, abc)
    if (censoring) {
      censored[censored %in% T] <- NA
      censored[censored %in% F] <- rexpit(-4 + 0.01 * age * abc)[censored %in% F]
      temp.df <- data.frame(temp.df, C = BinaryToCensoring(is.censored = censored))
      index <- index & censored %in% F
    }
    if (repeated.y || t == num.time.points) {
      prev.y <- y
      y <- rexpit(-3 + 0.003 * abc * age - 0.8 * abc * sex + 0.001 * (ldl - 120)^2)
      y[censored & !prev.y] <- NA
      y[prev.y == 1] <- 1
      temp.df <- data.frame(temp.df, Y = y)
      index <- index & y %in% 0
    }
    names(temp.df) <- paste0(names(temp.df), t)
    df <- data.frame(df, temp.df)
    
    prev.ldl <- ldl
    prev.abc <- abc
    prev.renal <- renal
    prev.y <- y
  }
  return(df)
}
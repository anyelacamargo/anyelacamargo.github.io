mh <- function(x0, f, dprop, rprop, N, B) {
  
  
  x <- matrix(NA, N + B, length(x0))
  fx <- rep(NA, N + B)
  x[1,] <- x0
  fx[1] <- f(x0)
  
  ct <- 0
  for(i in 2:(N + B)) {
    u <- rprop(x[i-1,])
    
    fu <- f(u)
    r <- log(fu) + log(dprop(x[i-1,], u)) - log(fx[i-1]) - log(dprop(u, x[i-1,]))
    R <- min(exp(r), 1)
    
    if(runif(1) <= R) {
      ct <- ct + 1
      x[i,] <- u
      fx[i] <- fu
    } else {
      x[i,] <- x[i-1,]
      fx[i] <- fx[i-1]
    }
    
    
  }
  
  return(list(x=x[-(1:B),], fx=fx[-(1:B)], rate=ct / (N + B)))
  
}


dpost <- function(theta) {
  
  a <- theta[1]
  b <- theta[2]
  p <- 1 - 1 / (1 + exp(a + b * x))
  lik <- exp(sum(dbinom(y, size=1, prob=p, log=TRUE)))
  dprior <- exp(a) * exp(-exp(a) / b.mme) / b.mme
  return(lik * dprior)
}


dprop <- function(theta, theta0) {
  
  a <- theta[1]
  b <- theta[2]
  pr1 <- exp(a) * exp(-exp(a) / b.mme) / b.mme
  pr2 <- dnorm(b, b.mle, sqrt(var.b.mle))
  return(pr1 * pr2)
}

rprop <- function(theta0) {
  a <- log(rexp(1, 1 / b.mme))
  b <- rnorm(1, b.mle, sqrt(var.b.mle))
  return(c(a, b))
}

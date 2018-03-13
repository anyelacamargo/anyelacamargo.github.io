

# Run Metropolis-Hastings
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

# Data


y <- c(1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0)
x <- c(53, 57, 58, 63, 66, 67, 67, 67, 68, 69, 70, 70, 70, 70, 72, 73, 75, 75,
       76, 76, 78, 79, 81)

# Preliminary output from ML estimation
logreg.out <- glm(y ~ x, family=binomial(logit))
summary(logreg.out)
pdata <- data.frame(x=seq(min(x), max(x),len=length(x)))
pdata$y = predict(logreg.out, newdata=pdata, type="response")
plot(x, y, pch=16,col="gray",
     ylab='Probability of O-rings Failing', xlab ='Temp')
lines(y ~ x, pdata, col="green4", lwd=2)

a.mle <- as.numeric(logreg.out$coefficients[1])
b.mle <- as.numeric(logreg.out$coefficients[2])
var.a.mle <- summary(logreg.out)$cov.scaled[1, 1]
var.b.mle <- summary(logreg.out)$cov.scaled[2, 2]
b.mme <- exp(a.mle + 0.577216)

# Posterior distribution

N <- 10000
B <- 1000
x0 <- c(a.mle, b.mle)
mh.out <- mh(x0, dpost, dprop, rprop, N, B)
alpha.mh <- mh.out$x[,1]
beta.mh <- mh.out$x[,2]

hist(alpha.mh, freq=FALSE, col="gray", border="white", xlab=expression(alpha))
plot(alpha.mh, type="l", col="gray", xlab="Iteration", ylab=expression(alpha))
lines(1:N, cumsum(alpha.mh) / (1:N))
hist(beta.mh, freq=FALSE, col="gray", border="white", xlab=expression(beta))
plot(beta.mh, type="l", col="gray", xlab="Iteration", ylab=expression(beta))
lines(1:N, cumsum(beta.mh)/(1:N))
p65 <- 1 - 1 / (1 + exp(alpha.mh + beta.mh * 65))
p45 <- 1 - 1 / (1 + exp(alpha.mh + beta.mh * 45))
hist(p65, freq=FALSE, col="gray", border="white", xlab="p(65)", main="")
hist(p45, freq=FALSE, col="gray", border="white", xlab="p(45)", main="")




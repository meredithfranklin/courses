## Newton's method to find the value of theta maximizing the likelihood (log likelihood)
## function. It does this by finding the root of the derivative of the log likelihood function.
## Arguments of newtonMethod function:
##  x - the sample data
##  d - first derivative of the likelihood function being maximized
##  d2 - second derivative of the likelihood function being maximized
##  theta0 - initial guess for theta
##  tol - tolerance level for convergence
##  maxiter - the maximum number of iterations for convergence
##
##  Newton's method approximates the root of a function, f, by the following
##  See slides for derivation from Taylor's series explansion
##     x_{n+1} = x_{n} - f''(x_{n}) / f'(x_{n})

library(MASS) # for generalized inverse ginv

newtonsMethod <- function(d, d2, theta0, tol=1e-08, maxiter= 1000) {
  theta <- theta0
  t<-0
  
  repeat{
    t<-t+1
    theta.hat <- theta - ginv(d2(theta)) %*% d(theta)
    if(mean(abs(theta.hat - theta)) < tol | t >= maxiter) {
      if(t >= maxiter) warning("Maximum number of iterations reached!")
      break
    }
    cat("Iteration:", t, ", theta.hat =", theta, "\n")
    theta<-theta.hat
  }  
  cat("Iteration:", t, ", theta.hat =", theta, "\n")
  out<-list(solution=theta.hat, n.iter=t, value=abs(theta.hat - theta))
  return(out)
}



## Gamma MLE - two parameter distribution
# Random sample of data with true alpha=5 (shape parameter), beta=2 (scale parameter)
# X<- rgamma(40, shape=5, scale = 2)
X <- c( 8.31,  6.54,  7.37, 10.87,  3.94, 10.12,  6.52, 13.41,  9.55,  6.18, 13.27,  9.33, 10.96,  4.92,
        11.88,  5.75, 17.75, 13.48,  6.93, 20.97,  7.40, 7.23,  4.45,  7.91, 16.06,  8.43, 10.55,  5.75,
        5.61, 10.74,  2.85,  6.74, 13.90,  6.94, 12.23,  3.58, 10.77, 20.19,  9.42, 10.07)

dl <- function(theta) {
  
  alpha <- theta[1]
  beta <- theta[2]
  n <- length(X)
  dl1 <- -n * log(beta) - n * digamma(alpha) + sum(log(X))
  dl2 <- -n * alpha / beta + n * mean(X) / beta**2
  return(c(dl1, dl2))
  
}

ddl <- function(theta) {
  
  alpha <- theta[1]
  beta <- theta[2]
  n <- length(X)
  ddl11 <- -n * trigamma(alpha)
  ddl12 <- -n / beta
  ddl22 <- -n * (2 * mean(X) / beta**3 - alpha / beta**2)
  return(matrix(c(ddl11, ddl12, ddl12, ddl22), 2, 2, byrow=TRUE))
  
}

newtonsMethod(dl, ddl, c(4.1,1.5))




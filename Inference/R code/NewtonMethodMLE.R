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
##     x_{n+1} = x_{n} - f(x_{n}) / f'(x_{n})

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

# Exponential MLE
# Random sample with theta=1
X<- rexp(40, rate = 1)

# First derivative
dl <- function(theta){
  -(length(X)/theta) + (sum(X)/theta**2)
}
# Second derivative
ddl <- function(theta){
  (length(X)/theta**2) - (2*sum(X)/theta**3)
}

# Apply Newton's Method
newtonsMethod(dl,ddl,0.5)

# Compare to sample mean
mean(X)

# Logistic MLE
# Random sample of data with true theta=5
X<- rlogis(40, location = 5, scale = 1)
#X <- c(6.40, 5.40, 4.30, 3.70, 4.67, 4.20, 5.74, 6.47, 6.19, 1.60, 5.30, 4.97, 5.29, 4.79,
 #      3.91, 4.78, 5.54, 4.90, 6.38, 6.13, 3.40, 4.69, 5.05, 5.46, 1.82, 2.24, 3.83, 4.46,
  #     7.66, 5.28, 2.69, 3.27, 5.37, 4.17, 6.42, 4.69, 8.43, 2.97, 8.26, 6.03)

dl <- function(theta) {
  length(X) - 2 * sum( exp(-(X - theta)) / (1 + exp(-(X - theta))))
}
ddl <- function(theta){
  -2 * sum( exp(-(X - theta)) / (1 + exp(-(X - theta)))**2 )
}

newtonsMethod(dl,ddl,mean(X))

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




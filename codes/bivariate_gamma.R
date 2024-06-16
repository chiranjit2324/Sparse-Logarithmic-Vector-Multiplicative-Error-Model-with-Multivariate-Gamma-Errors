
set.seed(12345)
library(mpoly)

#########################################
# Simulation from Bivariate Gamma:
#########################################
true_alpha = 2; true_beta = 4

yt1 = rgamma(n=100,shape=true_alpha,rate = true_beta)
yt2 = rgamma(n=100,shape=true_alpha,rate = true_beta)
yt3 = rgamma(n=100,shape=true_alpha,rate = true_beta)

xt1 = yt1 + yt2
xt2 = yt1 + yt3

# Data and dimensions:
X = data.frame(xt1,xt2)
d = ncol(X)

###################################################################################
# Function to calculate negative log-likelihood function (Using Numerical Method):
###################################################################################
neg_log_likelihood_numerical1 = function(param){
  alpha_par = param[1]; beta_par = param[2]
  integration_vector = vector()
  for(t in 1:nrow(X)){
    for(i in 1:d){
      integrand <- function(z){
        (X[t,i] - z)^(alpha_par - 1) * exp(-beta_par*(X[t,i] - z)) * z^(alpha_par - 1) * exp(-beta_par*z) * beta_par^(alpha_par)/gamma(alpha_par)
      }
      b_t <- min(X[t,])
    }
    integration_vector <- c(integration_vector,log(integrate(integrand,lower = 0,upper = b_t)$value))
  }
  N_t = nrow(X)
  term1 = N_t*(alpha_par * log(beta_par) - lgamma(alpha_par)) 
  term2 = sum(integration_vector)
  neg_log_lik = -sum(term1,term2)
  return(neg_log_lik)
}

###########################
# Calculate MLE:
###########################

# Bounds and max iterations for the parameters:
lb <- rep(1,2)
ub <- rep(10,2)
eps = 1e-5

# Initial value:
initial_value = rep(1.5,2)
neg_log_likelihood_numerical1(param = initial_value)

# Optimization:
library(nloptr)
bobyqa(x0 = initial_value,fn =  neg_log_likelihood_numerical1,
       lower = lb,upper = ub, control = list(maxeval = 3000))

library(Rsolnp)
solnp(pars = initial_value
      ,fun = neg_log_likelihood_numerical1,
      LB = lb,UB = ub, control = list(tol=eps,trace=1))


#########################################
# Simulation from Bivariate Gamma:
#########################################
true_alpha = 2; true_beta = 4

yt1 = rgamma(n=100,shape=true_alpha,rate = true_beta)
yt2 = rgamma(n=100,shape=true_alpha,rate = true_beta)
yt3 = rgamma(n=100,shape=true_alpha,rate = true_beta)

xt1 = yt1 + yt2
xt2 = yt1 + yt3

# Data and dimensions:
X = data.frame(xt1,xt2)
d = ncol(X)

log_factorial <- function(x) lgamma(x + 1)

neg_log_likelihood = function(param){
  #v1 = param[1];v2 = param[2];v3 = param[3]
  v1=v2=v3=param[1]
  beta1 = param[2]
  xt1 = X[,1]; xt2 = X[,2]
  sum1 = sum(log(dgamma(x = xt1,shape = v1+v2,rate = beta1)))  
  sum2 = sum(log(dgamma(x = xt2,shape = v1+v3,rate = beta1))) 
  
  new_func = function(r,v1,v2,v3,x1,x2) {
    factorial(r) * (gamma(v1)*gamma(v1+v2)*gamma(v1+v3))/(gamma(v1+r)*gamma(v1+v2+r)*gamma(v1+v3+r)) * as.function(laguerre(degree = r,alpha = v1+v2-1))(x1) * as.function(laguerre(degree = r,alpha = v1+v3-1))(x2)
  }
  
  N_t = nrow(X)
  ans_vec = rep(NA,N_t)
  for(i in 1:N_t){
    x1 = xt1[i]
    x2 = xt2[i]
    r <- 1
    tol <- 1e-4
    ans <- 0
    while (TRUE) {
      next_term <- new_func(r,v1,v2,v3,x1,x2)
      ans <- ans + next_term
      if (abs(next_term)<tol) break
      r <- r+1
      if(r>30) break
    }
    ans_vec[i] = log(1+ans)
  }
  
  sum3 = sum(ans_vec)  
  neg_log_lik = -(sum1+sum2+sum3)
  return(neg_log_lik)
}

###########################
# Calculate MLE:
###########################

# Bounds and max iterations for the parameters:
lb <- rep(1,2)
ub <- rep(10,2)
eps = 1e-5

# Initial value:
initial_value = rep(1.5,2)
neg_log_likelihood(param = initial_value)

# Optimization:
library(nloptr)
system.time({all_mle = bobyqa(x0 = initial_value,fn =  neg_log_likelihood,
       lower = lb,upper = ub, control = list(maxeval = 3000))})

library(nloptr)
system.time({all_mle = bobyqa(x0 = c(1,2),fn = neg_log_likelihood
                              ,lower = lb,upper = ub)})
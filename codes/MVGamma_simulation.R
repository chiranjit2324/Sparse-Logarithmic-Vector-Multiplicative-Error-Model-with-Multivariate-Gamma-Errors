
set.seed(12345)

#########################################
# Simulation from multivariate Gamma:
#########################################
simulate_mvgamma <- function(N,d,beta,alpha){
  lambda_vec = rep(alpha,d)
  theta_vec = rep(beta,d)  
  z <- rgamma(n=N,shape = alpha,rate = beta)
  v_mat <- matrix(NA,nrow=N,ncol = d)
  x_mat <- matrix(NA,nrow=N,ncol = d)
  for(i in 1:d){
    v_mat[,i] <- rgamma(n=N,shape =lambda_vec[i] ,rate =theta_vec[i])
    x_mat[,i] <- v_mat[,i] + z
  }
  return(x_mat)
}

X = simulate_mvgamma(N = 500,d = 3,beta=10,alpha=2)
d=ncol(X)

##########################################################
# Function to calculate negative log-likelihood function 
##########################################################
neg_log_likelihood_numerical1 = function(param){
  beta_par = param[1]; alpha_par = param[2]
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
eps = 1e-3

# Initial value:
initial_value = rep(1.5,2)
neg_log_likelihood_numerical1(param = initial_value)

# Optimization:
library(nloptr)
bobyqa(x0 = initial_value,fn =  neg_log_likelihood_numerical1,
       lower = lb,upper = ub, control = list(maxeval = 3000))

# library(Rsolnp)
# solnp(pars = initial_value
#       ,fun = neg_log_likelihood_numerical1,
#       LB = lb,UB = ub, control = list(tol=eps,trace=1))




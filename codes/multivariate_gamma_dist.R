
#############################################
# Function to calculate the pdf of mvgamma:
#############################################
mvgamma_pdf <- function(x_vec,alpha,beta,lambda_vec,theta_vec){
  m <- length(x_vec)
  
  # Constant:
  const <- (beta^alpha)/gamma(alpha)
  B <- min(x_vec)
  
  # Integrand function:
  integrand <- function(z){
    val <- (x - z)^(lambda - 1)* exp(-theta_vec[i]*(x - z)) * z^(alpha - 1) * exp(-beta*z)
  }
  
  prod_val <- rep(NA,m)
  for(i in 1:m){
    x <- x_vec[i] 
    lambda <- lambda_vec[i]
    integrand_val <- integrate(integrand,lower = 0,upper = B)
    prod_val[i] <- theta_vec[i]^(lambda_vec[i])/gamma(lambda_vec[i]) * integrand_val$value
  }
  
  # Final value of pdf:
  pdf_val <- prod(prod_val)*const
  return(pdf_val)
}

# #################
# # Parameters:
# #################
# m <- 3 #dimension
# alpha <- 2
# beta <- 1
# lambda_vec <- c(1,1,1)
# theta_vec <- c(1,1,1)
# 
# x_vec <- c(1,1,1)
# 
# # Evaluation of pdf:
# mvgamma_pdf(x_vec,alpha,beta,lambda_vec,theta_vec)

#########################################################################
# Simulate N observations from a m-dimensional multivariate Gamma vector
# of observations (Tsionas, 2003):
#################################################################

simulate_mvgamma <- function(N,m,alpha,beta,lambda_vec,theta_vec){
  z <- rgamma(n=N,shape = alpha,scale = 1/beta)
  v_mat <- matrix(NA,nrow=N,ncol = m)
  x_mat <- matrix(NA,nrow=N,ncol = m)
  for(i in 1:m){
    v_mat[,i] <- rgamma(n=N,shape =lambda_vec[i] ,scale =1/theta_vec[i])
    x_mat[,i] <- v_mat[,i] + z
  }
  return(x_mat)
}
  
# N <- 500
# alpha <- 1
# beta <- 2
# m <- 2
# lambda_vec <- rep(1,m)
# theta_vec <- rep(2,m)
# 
# x_mat <- simulate_mvgamma(N,m,alpha,beta,lambda_vec,theta_vec)
# 
# #######################
# # Correctness:
# #######################
# # Theoretical Mean:
# lambda_vec/theta_vec + alpha/beta
# 
# # Sample Mean:
# apply(x_mat,2,mean)
# 
# # Theoretical Covariance:
# cov_mat <- matrix(NA,nrow=m,ncol=m)
# for(i in 1:nrow(cov_mat)){
#   for(j in 1:ncol(cov_mat)){
#     if(i==j){
#       cov_mat[i,i] <- lambda_vec[i]/theta_vec[i]^2 + alpha/beta^2
#     }else{
#       cov_mat[i,j] <- alpha/beta^2
#     }
#   }
# }
# 
# cov_mat_theoretical <- cov_mat
# 
# #Sample Covariance:
# cov(x_mat)
# 





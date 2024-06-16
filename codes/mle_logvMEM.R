
# Load libraries:
library(ACDm)
library(nloptr)
library(Rsolnp)
library(dplyr)
library(doSNOW)

# Source codes:
source("codes/generate_autoregressive_matrices.R")
source("codes/sim_logvMEM.R")
source("codes/multivariate_gamma_dist.R")

##################################################
# Custom functions: Generate lags of a matrix:
##################################################
generate_lag_matrix <- function(x_t,t,lag_order){
  mat_list <- list()
  for(i in 1:lag_order){
    mat_list <- append(mat_list,list(x_t[t-i,]))
  }
  return(mat_list)
}

#generate_lag_matrix(x_t=X,t=10,lag_order=1)

###############################################
# Custom functions: Sum product of two lists:
###############################################

list.multiply <- function(list1,list2){
  a <- 0
  if(length(list1) == length(list2)){
    for(i in 1:length(list1)){
      m1 <- t(list1[[i]]) #matrix of dim d*d
      m2 <- as.matrix(unlist(list2[[i]]),ncol=1) # vector of dim d
      # a <- a+crossprod(list1[[i]],list2[[i]])
      a <- a + crossprod(m1,m2)
      
    }
  }
  return(a)
}

########################################
# Log-likelihood function for logvMEM:
########################################
log_lik_logvMEM <- function(X,A_matrix,beta_par){
  all_param <- unlist(lapply(1:ncol(A_matrix),function(col.index) A_matrix[,col.index]))
  N <- nrow(X)
  d <- ncol(X)
  
  # logvMEM parameters:
  logvMEM_par <- all_param
  A_matrices <- logvMEM_par[(1):(p*d^2)]
  # if(q==0){
  #   B_matrices <- NULL
  # }else{
  #   B_matrices <- logvMEM_par[(p*d^2+1):length(logvMEM_par)]
  # }
  
  A_matrix_list <- list()
  start_id <- seq(1,p*d^2,by=d^2)
  end_id <- seq(d^2,p*d^2,by=d^2)
  A_matrix_list <- lapply(1:p,function(i) matrix(A_matrices[start_id[i]:end_id[i]],nrow = d,ncol=d))
  
  # if(q==0){
  #   B_matrix_list <- NULL
  # }else{
  #   B_matrix_list <- list()
  #   start_id <- seq(1,q*d^2,by=d^2)
  #   end_id=seq(d^2,q*d^2,by=d^2)
  #   for(i in 1:q){
  #     B_matrix_list[[i]] <- matrix(B_matrices[start_id[i]:end_id[i]],nrow = d,ncol=d)
  #   }
  #
  # }
  
  # Generate mu_t:
  log_mu_t <- matrix(NA,nrow = N,ncol=d)
  log_mu_t[1:p,] <- rep(0.1,p) #initial values
  
  # Writing omega in terms of the elements of matrix A and data
  mu_vec <- apply(X,2,mean)
  variance_vec <- diag(cov(X))
  M <- rep(digamma(beta_par),d) - log(beta_par)
  omega <- crossprod(t((diag(x = 1,nrow = d,ncol = d) - Reduce(f = "+",x = A_matrix_list)))
                     , matrix((log(mu_vec) - variance_vec/(2*mu_vec^2)))) - M
  
  ################
  log_x_t <- log(X)
  
  log_mu_t[(p+1):N,] <- do.call(rbind,lapply((p+1):N, function(t) t(omega +
                                                                      list.multiply(list1 = A_matrix_list
                                                                                    ,list2 = generate_lag_matrix(x_t = log_x_t,t=t,lag_order = p)))))
  mu_t <- exp(log_mu_t)
  
  # Likelihood calculation:
  loglik1 <- -sum(log_mu_t[-c(1:p),])
  loglik2 <- ((beta_par/2)*log(beta_par) - log(gamma(beta_par/2))) * ((N - p - 2)*d + 1)
  
  a <- vector()
  for(i in 1:d){
    for(t in (p+1):N){
      x <- X[t,i]
      mu <- mu_t[t,i]
      
      integrand <- function(z){
        (x/mu - z)^(beta_par/2 - 1) * exp(-beta_par*(x/mu - z)) * beta_par^(beta_par/2)/gamma(beta_par/2) * z^(beta_par/2 - 1) * exp(-beta_par*z)
      }
      
      #b_t <- min(X[t,]/mu_t[t,])
      b_t <- x/mu
      #b_t <- 100
      a <- c(a,log(integrate(integrand,lower = 0,upper = b_t)$value))
    }
  }
  
  loglik3 <- sum(a)
  
  #loglik3 <- sum(data_matrix,na.rm = T)
  loglik <- sum(loglik1,loglik2,loglik3)
  return(loglik)
}

log_lik_logvMEM_ith_component <- function(X,A_matrix_i,i.index,beta_par){
  N <- nrow(X)
  d <- ncol(X)
  all_param_i <- A_matrix_i
  # # mvgamma parameters (fixed and known):
  # mvgamma_par <- c(alpha=beta_par/2,beta=beta_par,lambda_vec=rep(beta_par/2,d),theta_vec=rep(beta_par,d))
  # alpha <- mvgamma_par[1]
  # beta <- mvgamma_par[2]
  # lambda_vec <- mvgamma_par[3:(2+d)]
  # theta_vec <- mvgamma_par[(2+d+1):(2+2*d)]
  
  # logvMEM parameters:
  logvMEM_par <- all_param_i
  A_matrices <- logvMEM_par[1:(p*d)]
  # omega <- logvMEM_par[1]
  # if(p == 0){
  #   A_matrices <- NULL
  # }else{
  #   A_matrices <- logvMEM_par[1:(p*d)]
  # }
  # if(q==0){
  #   B_matrices <- NULL
  # }else{
  #   B_matrices <- logvMEM_par[(p*d+1):length(logvMEM_par)]
  # }
  
  # Writing omega in terms of the elements of matrix A and data
  start_id <- seq(1,length(A_matrices),by=d)
  end_id <- seq(d,length(A_matrices),by=d)
  A_matrices_list <- list()
  for(i in 1:length(start_id)){
    A_matrices_list[[i]] <- A_matrices[start_id[i]:end_id[i]]
  }
  mu_vec <- colMeans(X)
  variance_vec <- diag(cov(X))
  M <- rep(digamma(beta_par)- log(beta_par),d)[i.index] 
  omega <- tcrossprod((diag(x = 1,nrow = d)[i.index,] - Reduce(f = "+",x = A_matrices_list)) ,  t(matrix((log(mu_vec) - variance_vec/(2*mu_vec^2))))) - M
  
  log_x_t <- log(X)
  log_mu_i <- rep(0,N)
  for(t in (p+1):N){
    log_mu_i[t] <- omega + crossprod(A_matrices , matrix(unlist(generate_lag_matrix(x_t = log_x_t,t=t,lag_order = p))))
  }
  
  mu_i <- exp(log_mu_i)
  
  # Likelihood calculation based on ith component:
  loglik1 <- -sum(log_mu_i[-c(1:p)])
  loglik2 <- (((beta_par/2)*log(beta_par) - log(gamma(beta_par/2))) * ((N - p - 2)*d + 1))/d
  
  a <- vector()
#  for(i in 1:d){
    for(t in (p+1):N){
      x <- X[t,i.index]
      mu <- mu_i[t]
      
      integrand <- function(z){
        (x/mu - z)^(beta_par/2 - 1) * exp(-beta_par*(x/mu - z)) * beta_par^(beta_par/2)/gamma(beta_par/2) * z^(beta_par/2 - 1) * exp(-beta_par*z)
      }
      
      b_t <- x/mu
      #b_t <- x/mu
      #b_t <- 100
      a <- c(a,log(integrate(integrand,lower = 0,upper = b_t)$value))
    }
#  }
  
  loglik3 <- sum(a)
  loglik <- sum(loglik1,loglik2,loglik3)
  return(loglik)
}

#########################################
# Function to calculate initial values:
#########################################
init.val <- function(X,maxlag,n_beta){
  p <- maxlag
  d <- ncol(X)
  
  # All initial values:
  init_val_list <- list()
  par_matrix <- matrix(NA,nrow=d,ncol=p)
  par_matrix_list <- list()
  beta_i <- c()
  
  for(i in 1:d){
    lacd_fit <- acdFit(durations = X[,i], model = "LACD1",
                       dist = "gengamma", order = c(p,1),forceErrExpec = TRUE
                       ,startPara = c(0.1,rep(0.1,p),0,3,1),fixedParamPos = c(FALSE,rep(FALSE,p),TRUE,FALSE,TRUE))
    
    beta_i[i] <- lacd_fit$parameterInference["kappa",1]
    for(lag_id in 1:p){
      par_matrix[i,lag_id] <- lacd_fit$parameterInference[paste0("alpha",lag_id),1]
      par_matrix_list[[lag_id]] <- diag(par_matrix[,lag_id])
    }
  }
  
  beta_vec <- seq(ifelse(min(beta_i)-2>0,min(beta_i)-2,1),ifelse(max(beta_i)+2>0,max(beta_i+2),10)
                  ,length.out=n_beta)
  beta_vec_new <- seq(ifelse(min(beta_i)>0,min(beta_i),1),ifelse(max(beta_i)>0,max(beta_i),10)
                      ,length.out=n_beta)
  
  return(list(init_A_matrix = do.call(cbind,par_matrix_list),beta_list = beta_vec
              ,beta_list_max_min = beta_vec_new))
}

#################################
# Profile Likelihood estimation:
#################################

mle_ith_comp <- function(X,i.index,maxlag,init_A_matrix,beta_par,max_iter_bobyqa=3000,solver,eps=1e-3){
  p <- maxlag
  d <- ncol(X)
  
  # Bounds and max iterations for the parameters:
  lb <- rep(-0.5,p*d)
  ub <- rep(1,p*d)
  
  if(solver == "bobyqa"){  
    # Using bobyqa from nloptr:
    all_param_list <- tryCatch({bobyqa(x0 = init_A_matrix[i.index,]
                                       ,fn = function(param_i) -log_lik_logvMEM_ith_component(X=X,A_matrix_i = param_i
                                                                                              ,i.index = i.index,beta_par = beta_par),
                                       lower = lb,upper = ub, control = list(maxeval = max_iter_bobyqa,xtol_rel=eps))},
                               error=function(e) {
                                 return(NA)
                               })
    # Parameter estimates for each component:
    A_i_estimate <- all_param_list$par
    # Convergence for each component:
    conv_i <- all_param_list$convergence
    # Iteration for each component:
    iter_i <- all_param_list$iter
    # Log-likelihood for each component:
    loglik_i <- log_lik_logvMEM_ith_component(X = X,A_matrix_i = A_i_estimate
                                              ,beta_par = beta_par,i.index = i.index)
    
  }
  
  if(solver == "solnp"){
    # Using solnp from package Rsolnp:
    all_param_list <- tryCatch({solnp(pars = init_A_matrix[i.index,]
                                      ,fun = function(param_i) -log_lik_logvMEM_ith_component(X=X,A_matrix_i = param_i
                                                                                              ,i.index = i.index,beta_par = beta_par),
                                      LB = lb,UB = ub, control = list(tol=eps,trace=1))},
                               error=function(e) {
                                 return(NA)
                               }) 
    
    # Parameter estimates for each component:
    A_i_estimate <- all_param_list$pars
    # Convergence for each component:
    conv_i <- all_param_list$convergence
    # Iteration for each component:
    iter_i <- all_param_list$nfuneval
    # Log-likelihood for each component:
    loglik_i <- log_lik_logvMEM_ith_component(X = X,A_matrix_i = A_i_estimate
                                              ,beta_par = beta_par,i.index = i.index)
    
  }
  return(list(A_i_estimate = A_i_estimate,conv_i = conv_i,iter_i = iter_i
              ,loglik_i = loglik_i))
}


#####################Testing####################################


# ################################
# # Simulate data from logvMEM:
# ################################
# d <- 3
# p <- 3
# q <- 0
# r <- max(p,q)
# L <- matrix(c(rep(2,d),rep(3,d),rep(3,d)),nrow=d,ncol=d,byrow = TRUE) # HLag- componentwise
# A_matrix_list <- autoregressive_matrix(d = d,p = p,L = L,max_eigen_abs = 0.95)$A_matrix_list
# 
# X <- sim_logvMEM(N = 1500,d = d,omega=matrix(rep(0.1,d),ncol=1),A=A_matrix_list
#                  ,B=NULL,beta_par=4,burn_in = 500)
# A_matrix <- do.call(cbind,A_matrix_list)
# 
# # Plotting the time series and autocorrelation:
# plot.ts(X)
# cor(X,method = "kendall")

##################################################################################################
# Example to test whether the sum of individual lpglik is equal to the loglik of the full model:
##################################################################################################
# for(beta_par in 1:8){
#   l1= log_lik_logvMEM(X,A_matrix,beta_par)
#   print(l1)
#   
#   l2 = sum(unlist(lapply(1:d,function(i.index) log_lik_logvMEM_ith_component(X,A_matrix_i=A_matrix[i.index,],i.index=i.index,beta_par))))
#   print(l2)
# }

# # Testing:
# all_init_values <- init.val(X,maxlag=p,n_beta=10)
# init_A_matrix <- all_init_values$init_A_matrix
# beta_grid <- all_init_values$beta_list

# # Testing:
# mle_i_output <- mle_ith_comp(X = X,i.index = 2,maxlag = p,init_A_matrix = all_init_values$init_A_matrix
#                              ,max_iter_bobyqa = 3000,solver = "bobyqa",eps = 1e-3,beta_par = 4)

# Profile likelihood estimation:
# beta_i_index <- expand.grid(i_index = 1:d,beta = beta_grid)
# beta_grid_df <- data.frame(beta=beta_grid,model_id=1:length(beta_grid))
# beta_i_index_df <- inner_join(beta_i_index,beta_grid_df)
# 
# library(doParallel)
# 
# # Set up parallel estimation:
# no_cores <- detectCores()
# cl <- makeCluster(no_cores - 1)
# registerDoParallel(cl)
# registerDoSNOW(cl)
# 
# # print out the progress for every iteration
# progress <- function(n) cat(sprintf("task %d is complete\n", n))
# opts <- list(progress=progress)
# 
# 
# X_train = X
# 
# mle_output_list <- foreach(i=1:nrow(beta_i_index_df), .packages=c("ACDm","nloptr","Rsolnp"),
#                            .export = ls(globalenv()),.verbose = TRUE,.options.snow = opts) %dopar% {
#                              tryCatch({mle_ith_comp(X=X_train,i.index = beta_i_index_df$i_index[i],maxlag=p
#                                                     ,init_A_matrix = init_A_matrix,beta_par=beta_i_index_df$beta[i]
#                                                     ,max_iter_bobyqa = 3000,solver="bobyqa",eps=1e-3)}
#                                       ,error=function(e){return(NA)})
#                            }
# 
# stopCluster(cl)
# 
# # Binding all outputs:
# all_mle_list <- list()
# for(i in 1:length(mle_output_list)){
#   mle_list <- mle_output_list[[i]]
#   if(is.na(mle_list[1])==T){
#     all_mle_list[[i]] <- NA
#   }else{
#     mle_df = cbind.data.frame(beta_i_index_df[i,],t(mle_list$A_i_estimate),convergence = mle_list$conv_i
#                               ,n_iterations = mle_list$iter_i,loglik_i = mle_list$loglik_i)
#     colnames(mle_df) <- c(colnames(beta_i_index_df),
#                           paste("A",1:(p*d),sep="_"),"convergence","n_iterations"
#                           ,"loglik_i")
#     all_mle_list[[i]] <- mle_df
#   }
# }
# 
# # Row binding all mle results as data frame:
# all_mle_df <- do.call(rbind,all_mle_list)
# 
# # Take only the complete cases:
# all_mle_df <- all_mle_df[complete.cases(all_mle_df),]
# 
# # Filter out the models with incomplete information:
# split_all_mle_df_bymodelid <- split.data.frame(all_mle_df,f = all_mle_df$model_id)
# rm_id <- which(lapply(split_all_mle_df_bymodelid,nrow)<d)
# if(length(rm_id)>0){
#   split_all_mle_df_bymodelid <- split_all_mle_df_bymodelid[-rm_id]
# }
# 
# # Adding full data log-lik in the data frame:
# for(i in 1:length(split_all_mle_df_bymodelid)){
#   split_all_mle_df_bymodelid[[i]]$loglik_full_data <- sum(split_all_mle_df_bymodelid[[i]]$loglik_i)
# }
# 
# # Calculating log-likelihood by model_id:
# mle_id <- which.max(unlist(lapply(1:length(split_all_mle_df_bymodelid),function(i) unique(split_all_mle_df_bymodelid[[i]]$loglik_full_data))))
# A_matrix_estimated_mle <- split_all_mle_df_bymodelid[[mle_id]][,paste("A",1:(p*d),sep="_")]
# beta_mle <- unique(split_all_mle_df_bymodelid[[mle_id]][,"beta"])

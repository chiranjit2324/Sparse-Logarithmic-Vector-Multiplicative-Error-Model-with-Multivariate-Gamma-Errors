

###################
# Source codes:
###################

# 1. MLE: 
source("codes/mle_logvMEM.R")

# 2. Fitted values logvMEM:
source("codes/fitted_values_logvMEM.R")

# 3. Forecast function logvMEM:
source("codes/n_step_forecast_function_logvMEM.R")

####################
# R packages:
####################
library(nloptr)
library(ACDm)
library(Rsolnp)

###########################
# Simulated Data:
###########################
# d <- 3
# p <- 3
# q <- 0
# r <- max(p,q)
# L <- matrix(c(rep(2,d),rep(3,d),rep(3,d)),nrow=d,ncol=d,byrow = TRUE) # HLag- componentwise
# A_matrix_list <- autoregressive_matrix(d = d,p = p,L = L,max_eigen_abs = 0.95)$A_matrix_list
# 
# X <- sim_logvMEM(N = 1500,d = d,omega=matrix(rep(0.1,d),ncol=1),A=A_matrix_list
#                  ,B=NULL,beta_par=4,burn_in = 500)
# plot.ts(X)
# A <- do.call(cbind,A_matrix_list)
# # Plotting the time series and autocorrelation:
# plot.ts(X)
# cor(X)
# par(mfrow=c(3,1))
# lapply(1:d,function(i) acf(X[,i]))
# 
# A_matrix= do.call(cbind,A_matrix_list)
# max(A_matrix)
# max(abs(unlist(lapply(A_matrix_list,function(x) eigen(x)$values))))

# #############################
# # Simulated data:
# #############################
# nsim <- 1
# d <- 5
# p <- 10
# q <- 0
# r <- max(p,q)
# L <- matrix(c(rep(3,d),rep(4,d),rep(4,d),rep(5,d),rep(5,d)),nrow=d,ncol=d,byrow = TRUE) # HLag-componentwise
# 
# X <- list()
# A_matrix_list <- list()
# for(i in 1:nsim){
#   A_matrix_list[[i]] <- autoregressive_matrix(d = d,p = p,L = L,max_eigen_abs = 0.96)$A_matrix_list
#   X[[i]] <- sim_logvMEM(N = 1500,d = d,omega=matrix(rep(0.1,d),ncol=1),A=A_matrix_list[[i]]
#                         ,B=NULL,burn_in = 500,beta_par = 4)
# }
# 
# X <- X[[1]]
# A_matrix_list <- A_matrix_list[[1]]
# A_matrix <- do.call(cbind,A_matrix_list)


###########################
# Define the functions:
###########################

# Negative log-likelihood
f_i <- function(i.index,A_matrix_i,beta_par,N){
  return((-1/N)*log_lik_logvMEM_ith_component(X = X,A_matrix_i = A_matrix_i
                                              ,i.index = i.index,beta_par = beta_par))
}

# L2 norm of a vector:
norm_vec <- function(x) sqrt(sum(x^2))

# Testing:
# f_i(i.index = 1,A_matrix_i = A_matrix[1,],beta_par = 4,N = 1500)

################################################
# MLE for ith index/component
################################################

# mle_ith_comp <- function(X,i.index,p,d){
#   # Generation of initial values:
#   #library(ACDm)
#   init_val_list <- list()
#   for(i in 1:d){
#     lacd_fit <- acdFit(durations = X[,i], model = "LACD1",
#                        dist = "gengamma", order = c(p,1), dailyRestart = 1,forceErrExpec = TRUE)
#     init_val_list[[i]] <- diag(lacd_fit$parameterInference$Coef[2:(2+p-1)])
#     
#   }
#   # MLE of logvMEM:
#   init_i <- unlist(lapply(init_val_list,function(x) x[i.index,])) #initial values for MLE
#   lb <- rep(-0.5,p*d)
#   ub <- rep(1,p*d)
#   
#   #optimizer
#   # library(Rsolnp)
#   optimization <- solnp(pars = init_i,
#                         fun = function(x) -log_lik_logvMEM_ith_component(A_matrix_i = x,i.index,X = X),
#                         LB = lb,
#                         UB = ub,control = list(trace=0))
#   mle_i <- optimization$pars
#   ret_list <- list(mle_i=mle_i)
#   return(ret_list)
# }

# #################
# # Testing:
# #################
# # Parameters:
# i.index <- 1
# eps <- 10^(-4)
# p <- p
# max_iter <- 5000
# penalty <- "group-lasso"
# Hlag_structure <- "componentwise"
# verbose <- TRUE
# lambda_penalty <- 0.2
# all_mle <- mle_ith_comp(X = X,i.index = i.index,p,d)
# init_i <- mle_i <- all_mle$mle_i
# h=2
# true_par <- unlist(lapply(A_matrix_list,function(x) x[i.index,]))
# threshold <- 1e-4


#############################
# Penalty functions:
#############################

# SCAD penalty function:
scad <- function(theta,gamma_par,lambda_par){
  gamma_lambda <- gamma_par*lambda_par
  if(theta <= lambda_par){
    val <- lambda_par*theta
  }
  if(lambda_par <= theta & theta <= gamma_lambda){
    val <- (gamma_par*lambda_par*theta - 0.5*(theta^2 + lambda_par^2))/(gamma_par - 1)
  }
  if(theta > gamma_lambda){
    val <- lambda_par^2*(gamma_par + 1)/2
  }
  return(val)
}

# Minimax Concave Penalty:
mcp <- function(theta,gamma_par,lambda_par){
  gamma_lambda <- gamma_par*lambda_par
  if(abs(theta)<=gamma_lambda){
    val <- lambda_par*abs(theta) - theta^2/(2*gamma_par)
  }else{
    val <- gamma_par*lambda_par^2/2
  }
  return(val)
}

# Penalty:
g_i <- function(i.index,A_matrix_i,lambda_penalty,penalty,Hlag_structure,mle_i_vec){
  
  ####################
  # Componentwise:
  ####################
  if(penalty == "group-lasso" & Hlag_structure == "componentwise"){
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    ret_value <- sum(unlist(lapply(group_list,function(x) norm_vec(A_mat[x]))))*lambda_penalty
  }
  if(penalty == "adaptive-group-lasso" & Hlag_structure == "componentwise"){
    gamma_par <- 1
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    A_mat_norm_groupwise <- unlist(lapply(group_list,function(x) norm_vec(A_mat[x])))
    mle_i_norm_groupwise <- unlist(lapply(group_list,function(x) norm_vec(mle_i_vec[x])))
    ret_value <- sum(mle_i_norm_groupwise^(-gamma_par)*A_mat_norm_groupwise)*lambda_penalty
  }
  if(penalty == "group-scad" & Hlag_structure == "componentwise"){
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    ret_value <- sum(unlist(lapply(group_list,function(x) scad(norm_vec(A_mat[x])
                                                               ,gamma_par = 3.7,lambda_par = lambda_penalty))))
  }
  if(penalty == "group-mcp" & Hlag_structure == "componentwise"){
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    ret_value <- sum(unlist(lapply(group_list,function(x) mcp(norm_vec(A_mat[x])
                                                              ,gamma_par = 3,lambda_par = lambda_penalty))))
  }
  
  #####################
  # Elementwise:
  #####################
  if(penalty == "group-lasso" & Hlag_structure == "elementwise"){
    A_mat_list <- lapply(seq(1,p*d,d),function(x) as.numeric(A_matrix_i[x:(x+d-1)]))
    # New formulation of A_matrix, rows correspond to lag_index and columns are column of A_mat
    new_A_mat <- do.call(rbind,A_mat_list)
    grouped_cal <- rep(NA,ncol(new_A_mat))
    for(i in 1:ncol(new_A_mat)){
      grouped_cal[i] <- sum(unlist(lapply(seq(1,p), function(x) norm_vec(new_A_mat[,i][x:p]))))
    }
    ret_value <- sum(grouped_cal)*lambda_penalty
  }
  if(penalty == "adaptive-group-lasso" & Hlag_structure == "elementwise"){
    gamma_par <- 1
    A_mat_list <- lapply(seq(1,p*d,d),function(x) as.numeric(A_matrix_i[x:(x+d-1)]))
    mle_i_list <- lapply(seq(1,p*d,d),function(x) as.numeric(mle_i_vec[x:(x+d-1)]))
    # New formulation of A_matrix, rows correspond to lag_index and columns are column of A_mat
    new_A_mat <- do.call(rbind,A_mat_list)
    new_mle_i_mat <- do.call(rbind,mle_i_list)
    grouped_cal <- rep(NA,ncol(new_A_mat))
    id <- which(new_mle_i_mat<0.005,arr.ind = T)
    new_mle_i_mat[id] <- 0.005
    for(i in 1:ncol(new_A_mat)){
      grouped_cal[i] <- sum(unlist(lapply(seq(1,p), function(x) norm_vec(new_A_mat[,i][x:p])*norm_vec(new_mle_i_mat[,i][x:p])^(-gamma_par))))
      #print(unlist(lapply(seq(1,p), function(x) as.matrix(new_A_mat[,i][x:p]))))
      #print(unlist(lapply(seq(1,p), function(x) as.matrix(new_mle_i_mat[,i][x:p]))))
    }
    ret_value <- sum(grouped_cal)*lambda_penalty
  }
  if(penalty == "group-scad" & Hlag_structure == "elementwise"){
    A_mat_list <- lapply(seq(1,p*d,d),function(x) as.numeric(A_matrix_i[x:(x+d-1)]))
    # New formulation of A_matrix, rows correspond to lag_index and columns are column of A_mat
    new_A_mat <- do.call(rbind,A_mat_list)
    grouped_cal <- rep(NA,ncol(new_A_mat))
    for(i in 1:ncol(new_A_mat)){
      grouped_cal[i] <- sum(unlist(lapply(seq(1,p), function(x) scad(norm_vec(new_A_mat[,i][x:p])
                                                                     ,gamma_par = 3.7,lambda_par = lambda_penalty))))
    }
    ret_value <- sum(grouped_cal)
  }
  if(penalty == "group-mcp" & Hlag_structure == "elementwise"){
    A_mat_list <- lapply(seq(1,p*d,d),function(x) as.numeric(A_matrix_i[x:(x+d-1)]))
    # New formulation of A_matrix, rows correspond to lag_index and columns are column of A_mat
    new_A_mat <- do.call(rbind,A_mat_list)
    grouped_cal <- rep(NA,ncol(new_A_mat))
    for(i in 1:ncol(new_A_mat)){
      grouped_cal[i] <- sum(unlist(lapply(seq(1,p), function(x) mcp(norm_vec(new_A_mat[,i][x:p])
                                                                    ,gamma_par = 3,lambda_par = lambda_penalty))))
    }
    ret_value <- sum(grouped_cal)
  }
  
  ###########################
  # Own-other:
  ###########################
  if(penalty == "group-lasso" & Hlag_structure == "own-other"){
    # Componentwise sum:
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    componentwise_sum <- sum(unlist(lapply(group_list,function(x) norm_vec(A_mat[x]))))
    # Own other sum:
    A_mat_list <- lapply(seq(1,p*d,d),function(x) A_matrix_i[x:(x+d-1)])
    new_l <- list()
    for(l in 1:(p-1)){
      new_l[[l]] <- norm_vec(c(as.numeric(A_mat_list[[l]][-i.index])
                                     ,as.numeric(sapply((l+1):p,function(x) A_mat_list[[x]]))))
    }
    own_other_sum <- sum(unlist(new_l))
    fun_val <- (componentwise_sum + own_other_sum) * lambda_penalty 
    ret_value <- fun_val
  }
  if(penalty == "adaptive-group-lasso" & Hlag_structure == "own-other"){
    # Componentwise sum:
    gamma_par <- 1
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    A_mat_norm_groupwise <- unlist(lapply(group_list,function(x) norm_vec(A_mat[x])))
    mle_i_norm_groupwise <- unlist(lapply(group_list,function(x) norm_vec(mle_i_vec[x])))
    componentwise_sum <- sum(mle_i_norm_groupwise^(-gamma_par)*A_mat_norm_groupwise)
    
    # Own other sum:
    A_mat_list <- lapply(seq(1,p*d,d),function(x) A_matrix_i[x:(x+d-1)])
    id <- which(mle_i_vec<0.005)
    mle_i_vec[id] <- 0.005
    
    mle_i_mat_list <- lapply(seq(1,p*d,d),function(x) mle_i_vec[x:(x+d-1)])
    new_l <- list()
    new_mle_l <- list()
    ratio_new_l_new_mle_l <- list()
    for(l in 1:(p-1)){
      new_l[[l]] <- norm_vec(c(as.numeric(A_mat_list[[l]][-i.index])
                                     ,as.numeric(sapply((l+1):p,function(x) A_mat_list[[x]]))))
      new_mle_l[[l]] <- norm_vec(c(as.numeric(mle_i_mat_list[[l]][-i.index])
                                         ,as.numeric(sapply((l+1):p,function(x) mle_i_mat_list[[x]]))))
      ratio_new_l_new_mle_l[[l]] <- new_mle_l[[l]]^(-gamma_par) * new_l[[l]]
    }
    own_other_sum <- sum(unlist(ratio_new_l_new_mle_l))
    fun_val <- (componentwise_sum + own_other_sum) * lambda_penalty 
    ret_value <- fun_val
  }
  if(penalty == "group-scad" & Hlag_structure == "own-other"){
    # Componentwise sum:
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    componentwise_sum <- sum(unlist(lapply(group_list,function(x) scad(norm_vec(A_mat[x])
                                                                       ,gamma_par = 3.7,lambda_par = lambda_penalty))))
    # Own other sum:
    A_mat_list <- lapply(seq(1,p*d,d),function(x) A_matrix_i[x:(x+d-1)])
    new_l <- list()
    for(l in 1:(p-1)){
      new_l[[l]] <- scad(norm_vec(c(as.numeric(A_mat_list[[l]][-i.index])
                                          ,as.numeric(sapply((l+1):p,function(x) A_mat_list[[x]]))))
                         ,gamma_par = 3.7,lambda_par = lambda_penalty)
    }
    own_other_sum <- sum(unlist(new_l))
    fun_val <- (componentwise_sum + own_other_sum)
    ret_value <- fun_val
  }
  if(penalty == "group-mcp" & Hlag_structure == "own-other"){
    # Componentwise sum:
    A_mat <- A_matrix_i
    group_list <- lapply(seq(1,d*p,d),function(x) x:(d*p))
    componentwise_sum <- sum(unlist(lapply(group_list,function(x) mcp(norm_vec(A_mat[x])
                                                                      ,gamma_par = 3,lambda_par = lambda_penalty))))
    # Own other sum:
    A_mat_list <- lapply(seq(1,p*d,d),function(x) A_matrix_i[x:(x+d-1)])
    new_l <- list()
    for(l in 1:(p-1)){
      new_l[[l]] <- mcp(norm_vec(c(as.numeric(A_mat_list[[l]][-i.index])
                                         ,as.numeric(sapply((l+1):p,function(x) A_mat_list[[x]]))))
                        ,gamma_par = 3,lambda_par = lambda_penalty)
    }
    own_other_sum <- sum(unlist(new_l))
    fun_val <- (componentwise_sum + own_other_sum)
    ret_value <- fun_val
  }
  return(ret_value)
}

# Function to be minimized:
F_i <- function(i.index,A_matrix_i,lambda_penalty,penalty,Hlag_structure,N,mle_i_vec,beta_par){
  f_val <- f_i(i.index,A_matrix_i,beta_par,N)
  g_val <- g_i(i.index,A_matrix_i,lambda_penalty,penalty,Hlag_structure,mle_i_vec)
  return(f_val + g_val)
}

# #####################################################################
# # Penalized likelihood estimation for fixed beta and fixed lambda: 
# # For d=5, p = 10, N = 1500, this takes about 136s for one iteration 
# #####################################################################
# 
# penalized_logvMEM_fix_beta_fix_lambda <- function(X_train,X_test,maxlag,max_iter,eps,penalty,Hlag_structure
#                                                    ,verbose,lambda_penalty,init_A_i_list,threshold,mle_A_matrix,beta_par,h_step){
#   X <- X_train
#   
#   # Max-lag:
#   p <- maxlag
#   # No. of observations:
#   N <- nrow(X)
#   
#   #############.
#   # Optimizer:
#   #############
#   library(nloptr)
# 
#   # Bounds and max iterations for the parameters:
#   lb <- rep(-0.5,p*d)
#   ub <- rep(1,p*d)
#   
#   A_list <- list()
#   neg_pen_loglik_vec <- c()
#   conv_vec <- c()
#   iteration_vec <- c()
#   
#   library(doParallel)
#   
#   # Set up parallel estimation:
#   no_cores <- detectCores()
#   cl <- makeCluster(no_cores - 1)
#   registerDoParallel(cl)
#   
#   A_matrix_i_list <- foreach(i.index=1:d, .packages=c("ACDm","nloptr"),.export = ls(globalenv()),.verbose = verbose) %dopar% {
#     tryCatch({bobyqa(x0 = unlist(init_A_i_list[[i.index]])
#                      ,fn = function(x) F_i(i.index = i.index, A_matrix_i = x
#                                            ,lambda_penalty =lambda_penalty,penalty=penalty,Hlag_structure = Hlag_structure
#                                            ,N=N,mle_i_vec=mle_A_matrix[i.index,]
#                                            ,beta_par = beta_par),
#                      lower = lb,upper = ub, control = list(maxeval = max_iter,xtol_rel=eps))},
#              error=function(e) {
#                return(NA)
#              })
#   }
#   
#   stopCluster(cl)
#   ##########################
#   # Components Metrics:  
#   ##########################
#   
#   # Parameter estimates for each component:
#   A_i_estimate_list <- lapply(1:d,function(i.index) A_matrix_i_list[[i.index]]$par)
#   # Convergence for each component:
#   conv_vec <- unlist(lapply(1:d,function(i.index) A_matrix_i_list[[i.index]]$convergence))
#   # Iteration for each component:
#   iter_vec <- unlist(lapply(1:d,function(i.index) A_matrix_i_list[[i.index]]$iter))
#   
#   # Penalized log likelihood for each component:
#   pen_loglik_vec <- unlist(lapply(1:d, function(i.index) -F_i(i.index = i.index, A_matrix_i = A_i_estimate_list[[i.index]]
#                                                                  ,lambda_penalty =lambda_penalty,penalty=penalty
#                                                                  ,Hlag_structure = Hlag_structure
#                                                                  ,N=N,mle_i_vec=mle_A_matrix[i.index,], beta_par = beta_par)))
#   loglik_i_vec <- unlist(lapply(1:d,function(i.index) log_lik_logvMEM_ith_component(X = X,A_matrix_i = A_i_estimate_list[[i.index]]
#                                                                                     ,beta_par = beta_par,i.index = i.index)))
#   
#   # AIC and BIC (Component):
#   df_i_vec <- unlist(lapply(1:d,function(i.index)length(which(abs(A_i_estimate_list[[i.index]])>threshold))))
#   aic_i_vec <- unlist(lapply(1:d,function(i.index) 2*df_i_vec[i.index] - 2*loglik_i_vec[i.index]))
#   bic_i_vec <- unlist(lapply(1:d,function(i.index) df_i_vec[i.index]*log(N) - 2*loglik_i_vec[i.index]))
#   
#   #######################################################
#   # Calculate Mean absolute deviation (Full data):
#   #######################################################
#   fitted_values_componentwise <- lapply(1:d, function(i.index) fitted_values_logvMEM_i(X = X,A_matrix_i = A_i_estimate_list[[i.index]]
#                                                                                        ,i.index = i.index,beta_par = beta_par,p = p,q=q))
# 
#   mad_i <- unlist(lapply(1:d,function(i.index) mean(abs(fitted_values_componentwise[[i.index]] - X[,i.index]))))
#   count_non_zero_coef_i=df_i_vec
#   
#   all_output <- list(A_list = A_i_estimate_list,beta_pars = beta_par,conv_vec=conv_vec,
#                      iter_vec = iter_vec,neg_pen_loglik_vec = neg_pen_loglik_vec,
#                      loglik_i_vec = loglik_i_vec,df_i_vec=df_i_vec,aic_i_vec=aic_i_vec
#                      ,bic_i_vec=bic_i_vec,mad_i)
#   
#   #####################
#   # Full data Metrics:
#   #####################
#   A_matrix_estimated_vectorized <- unlist(lapply(1:length(A_i_estimate_list),function(i) A_i_estimate_list[[i]]))
#   A_matrix_estimated <- do.call(rbind,A_i_estimate_list)
#   loglik_full_data <- log_lik_logvMEM(X = X,A_matrix = A_matrix_estimated,beta_par = beta_par)
#   
#   # AIC and BIC (Full data):
#   df_full_data <- length(which(abs(A_matrix_estimated_vectorized)>threshold))
#   aic_full_data <- 2*df_full_data - 2*loglik_full_data
#   bic_full_data <- df_full_data*log(N) - 2*loglik_full_data 
# 
#   mad_full_data <- mean(mad_i)
#   count_non_zero_coef_full_data=df_full_data
#   pen_loglik_full_data <- sum(pen_loglik_vec)
#   
#   # Omega 
#   mu_vec <- apply(X,2,mean)
#   
#   variance_vec <- diag(cov(X))
#   M <- rep(digamma(beta_par),d) - log(beta_par)
#   
#   # Convert A_matrix to list of A matrices:
#   A_matrix_list <- list()
#   start_id <- seq(1,p*d^2,by=d^2)
#   end_id <- seq(d^2,p*d^2,by=d^2)
#   A_matrix <- do.call(rbind,A_i_estimate_list)
#   for(i in 1:p){
#     A_matrix_list[[i]] <- matrix(A_matrix[start_id[i]:end_id[i]],nrow = d,ncol=d)
#   }
#   
#   omega <- crossprod(t((diag(x = 1,nrow = d,ncol = d) - Reduce(f = "+",x = A_matrix_list)))
#                      , matrix((log(mu_vec) - variance_vec/(2*mu_vec^2)))) - M
#   
#   # h-step forecast metrics:
#   forecast_metrics_logvMEM_full_data <- forecast_logvMEM_full_data(X_train=X_train,X_test=X_test,h_step=h_step
#                                          ,omega_estimated=omega,A_matrix_estimated=A_matrix_estimated,maxlag=maxlag)
#   
#   
#   all_output_df <- data.frame(component = 1:d,beta_par,lambda_penalty,omega,do.call(rbind,A_i_estimate_list)
#                               ,conv_vec,iter_vec,pen_loglik_vec
#                               ,loglik_i_vec,df_i_vec,aic_i_vec,bic_i_vec,penalty,Hlag_structure
#                               ,mad_i,count_non_zero_coef_i,loglik_full_data,df_full_data,aic_full_data,bic_full_data
#                               ,mad_full_data,count_non_zero_coef_full_data,pen_loglik_full_data
#                               ,as.numeric(forecast_metrics_logvMEM_full_data)[1]
#                               ,as.numeric(forecast_metrics_logvMEM_full_data)[2]
#                               ,as.numeric(forecast_metrics_logvMEM_full_data)[3])
# 
#   colnames(all_output_df) <- c("component","beta","lambda","omega",paste("A",1:(d*p),sep="_")
#                                ,"convergence_i", "n_iterations_i", "pen_loglik_i"
#                                ,"loglik_i", "df_i", "aic_i", "bic_i_vec","penalty"
#                                ,"HLag_structure", "MAD_i","count_non_zero_coef_i"
#                                ,"loglik_full_data","df_full_data"
#                                ,"aic_full_data","bic_full_data","MAD_full_data_insample"
#                                ,"count_non_zero_coef_full_data","pen_loglik_full_data"
#                                ,names(forecast_metrics_logvMEM_full_data))
#   
#   return(all_output_df)
# }

# #####################################################################
# # Penalized likelihood estimation for fixed beta and fixed lambda: 
# # For d=5, p = 10, N = 1500, this takes about 136s for one iteration 
# #####################################################################
penalized_logvMEM_fix_beta_fix_lambda <- function(X_train,X_test,i.index,maxlag,max_iter_bobyqa,eps,penalty,Hlag_structure
                                                  ,verbose,lambda_penalty,init_A_i_list,threshold,mle_A_matrix,beta_par,h_step
                                                  ,solver){
  X <- X_train
  # Max-lag:
  p <- maxlag
  # No. of observations:
  N <- nrow(X)
  
  #############################################
  # Parameters specification for optimization:
  #############################################
  library(nloptr)
  
  # Bounds and max iterations for the parameters:
  lb <- rep(-0.5,p*d)
  ub <- rep(1,p*d)
  
  if(solver == "bobyqa"){  
    # Using bobyqa from nloptr:
    A_matrix_i_list <- tryCatch({bobyqa(x0 = unlist(init_A_i_list[[i.index]])
                                        ,fn = function(x) F_i(i.index = i.index, A_matrix_i = x
                                                              ,lambda_penalty =lambda_penalty,penalty=penalty
                                                              ,Hlag_structure = Hlag_structure
                                                              ,N=N,mle_i_vec=mle_A_matrix[i.index,]
                                                              ,beta_par = beta_par),
                                        lower = lb,upper = ub, control = list(maxeval = max_iter_bobyqa,xtol_rel=eps))},
                                error=function(e) {
                                  return(NA)
                                })
    # Parameter estimates for each component:
    A_i_estimate <- A_matrix_i_list$par
    # Convergence for each component:
    conv_i <- A_matrix_i_list$convergence
    # Iteration for each component:
    iter_i <- A_matrix_i_list$iter
    # Penalized log likelihood for each component:
    pen_loglik_i <- -F_i(i.index = i.index, A_matrix_i = A_i_estimate
                         ,lambda_penalty =lambda_penalty,penalty=penalty
                         ,Hlag_structure = Hlag_structure
                         ,N=N,mle_i_vec=mle_A_matrix[i.index,], beta_par = beta_par)
    # Log-likelihood for each component:
    loglik_i <- log_lik_logvMEM_ith_component(X = X,A_matrix_i = A_i_estimate
                                              ,beta_par = beta_par,i.index = i.index)
    # AIC and BIC (Component):
    df_i <- length(which(abs(A_i_estimate)>threshold))
    aic_i <- 2*df_i - 2*loglik_i
    bic_i <- df_i*log(N) - 2*loglik_i
    # Writing omega in terms of the elements of matrix A and data
    start_id <- seq(1,length(A_i_estimate),by=d)
    end_id <- seq(d,length(A_i_estimate),by=d)
    A_matrices_list <- list()
    for(i in 1:length(start_id)){
      A_matrices_list[[i]] <- A_i_estimate[start_id[i]:end_id[i]]
    }
    mu_vec <- colMeans(X)
    variance_vec <- diag(cov(X))
    M <- rep(digamma(beta_par)- log(beta_par),d)[i.index] 
    omega_i <- (diag(x = 1,nrow = d)[i.index,] - Reduce(f = "+",x = A_matrices_list)) %*% matrix((log(mu_vec) - variance_vec/(2*mu_vec^2))) - M
    # Calculate the fitted values:
    fitted_values_componentwise <- fitted_values_logvMEM_i(X = X,A_matrix_i = A_i_estimate
                                                           ,i.index = i.index,beta_par = beta_par,p = p,q=q)
    # Calculate MAD:
    mad_i <- mean(abs(fitted_values_componentwise - X[,i.index]))
    # Count no. of nonzero coefficients:
    count_non_zero_coef_i=df_i
  }
  
  if(solver == "solnp"){
    # Using solnp from package Rsolnp:
    A_matrix_i_list <- tryCatch({solnp(pars = unlist(init_A_i_list[[i.index]])
                                       ,fun = function(x) F_i(i.index = i.index, A_matrix_i = x
                                                              ,lambda_penalty =lambda_penalty,penalty=penalty,Hlag_structure = Hlag_structure
                                                              ,N=N,mle_i_vec=mle_A_matrix[i.index,]
                                                              ,beta_par = beta_par),
                                       LB = lb,UB = ub, control = list(tol=eps,trace=1))},
                                error=function(e) {
                                  return(NA)
                                })  
    # Parameter estimates for each component:
    A_i_estimate <- A_matrix_i_list$pars
    # Convergence for each component:
    conv_i <- A_matrix_i_list$convergence
    # Iteration for each component:
    iter_i <- A_matrix_i_list$nfuneval
    # Penalized log likelihood for each component:
    pen_loglik_i <- -F_i(i.index = i.index, A_matrix_i = A_i_estimate
                         ,lambda_penalty =lambda_penalty,penalty=penalty
                         ,Hlag_structure = Hlag_structure
                         ,N=N,mle_i_vec=mle_A_matrix[i.index,], beta_par = beta_par)
    # Log-likelihood for each component:
    loglik_i <- log_lik_logvMEM_ith_component(X = X,A_matrix_i = A_i_estimate
                                              ,beta_par = beta_par,i.index = i.index)
    # AIC and BIC (Component):
    df_i <- length(which(abs(A_i_estimate)>threshold))
    aic_i <- 2*df_i - 2*loglik_i
    bic_i <- df_i*log(N) - 2*loglik_i
    # Writing omega in terms of the elements of matrix A and data
    start_id <- seq(1,length(A_i_estimate),by=d)
    end_id <- seq(d,length(A_i_estimate),by=d)
    A_matrices_list <- list()
    for(i in 1:length(start_id)){
      A_matrices_list[[i]] <- A_i_estimate[start_id[i]:end_id[i]]
    }
    mu_vec <- colMeans(X)
    variance_vec <- diag(cov(X))
    M <- rep(digamma(beta_par)- log(beta_par),d)[i.index] 
    omega_i <- (diag(x = 1,nrow = d)[i.index,] - Reduce(f = "+",x = A_matrices_list)) %*% matrix((log(mu_vec) - variance_vec/(2*mu_vec^2))) - M
    # Calculate the fitted values:
    fitted_values_componentwise <- fitted_values_logvMEM_i(X = X,A_matrix_i = A_i_estimate
                                                           ,i.index = i.index,beta_par = beta_par,p = p,q=q)
    # Calculate MAD:
    mad_i <- mean(abs(fitted_values_componentwise - X[,i.index]))
    # Count no. of nonzero coefficients:
    count_non_zero_coef_i=df_i
  }
  
  # List of all outputs:
  all_output <- list(A_list = A_i_estimate,beta_pars = beta_par,convergence=conv_i,
                     n_iter = iter_i, pen_loglik_i = pen_loglik_i,
                     loglik_i = loglik_i,df_i=df_i,aic_i=aic_i
                     ,bic_i=bic_i,omega_i = omega_i,mad_i =mad_i)
  
  # Output data frame:
  all_output_df <- data.frame(component = i.index,beta_par,lambda_penalty,omega_i,t(A_i_estimate)
                              ,conv_i,iter_i,pen_loglik_i
                              ,loglik_i,df_i,aic_i,bic_i,penalty,Hlag_structure
                              ,mad_i,count_non_zero_coef_i)
  # Assigning column names to output data frame:
  colnames(all_output_df) <- c("component","beta","lambda","omega_i",paste("A",1:(d*p),sep="_")
                               ,"convergence_i", "n_iterations_i", "pen_loglik_i"
                               ,"loglik_i", "df_i", "aic_i", "bic_i","penalty"
                               ,"HLag_structure", "MAD_i_insample","count_non_zero_coef_i")
  
  return(all_output_df)
}

# ##############
# # Testing:
# #############
# 
# # Data and parameters:
# X_train=X[1:1499,];X_test = X[1500,];h_step = 1;maxlag=p;max_iter_bobyqa=10;eps=10^(-5)
# penalty = "adaptive-group-lasso";Hlag_structure = "componentwise"
# verbose = TRUE;lambda_penalty = 0.01; threshold =1e-4
# all_initial_values <- init.val(X,maxlag=p,n_beta=1)
# init_A_matrix <- all_initial_values$init_A_matrix
# beta_grid <- all_initial_values$beta_list
# solver="bobyqa";eps=1e-3
# 
# # Set up parallel estimation:
# library(doParallel)
# no_cores <- detectCores()
# cl <- makeCluster(no_cores - 1)
# registerDoParallel(cl)
# registerDoSNOW(cl)
# 
# # print out the progress for every iteration
# progress <- function(n) cat(sprintf("task %d is complete\n", n))
# opts <- list(progress=progress)
# beta_i_index <- expand.grid(i_index = 1:d,beta = beta_grid)
# beta_grid_df <- data.frame(beta=beta_grid,model_id=1:length(beta_grid))
# beta_i_index_df <- inner_join(beta_i_index,beta_grid_df)
# 
# mle_output_list <- foreach(i=1:nrow(beta_i_index_df), .packages=c("ACDm","nloptr","Rsolnp"),
#                            .export = ls(globalenv()),.verbose = TRUE,.options.snow = opts) %dopar% {
#                              tryCatch({mle_ith_comp(X=X_train,i.index = beta_i_index_df$i_index[i],maxlag=p
#                                                     ,init_A_matrix = init_A_matrix,beta_par=beta_i_index_df$beta[i]
#                                                     ,max_iter_bobyqa = max_iter_bobyqa,solver=solver,eps=eps)}
#                                       ,error=function(e){return(NA)})
#                            }
# 
# stopCluster(cl)
# 
# 
# mle_A_matrix <- do.call(rbind,lapply(1:length(mle_output_list), function(x) mle_output_list[[x]]$A_i_estimate))
# 
# system.time({penalized_estimate_list <- penalized_logvMEM_fix_beta_fix_lambda(X_train = X_train,X_test = X_test,h_step = h_step,maxlag=p
#                                                         ,max_iter_bobyqa=2000,eps=10^(-3)
#                                                         ,penalty="group-lasso",Hlag_structure="elementwise"
#                                                         ,verbose = TRUE,lambda_penalty=0.01
#                                                         ,init_A_i_list=lapply(1:d,function(i.index) all_initial_values$init_A_matrix[i.index,])
#                                                         ,threshold=1e-4,mle_A_matrix=mle_A_matrix
#                                                         ,beta_par=4,solver = solver,i.index = 1)
# })
# 
# penalized_estimate_list


rm(list=ls())

# Load libraries
library(dplyr)
library(ACDm)
library(doParallel)
library(nloptr)
library(Rsolnp)
library(doSNOW)

# Start time:
start_time <- Sys.time()

# Source codes:
source("codes/sparse_logvMEM_estimation.R")
source("codes/sparsity_plot.R")

# Global parameters that control the optimization:
eps = 1e-3;threshold = 1e-3
max_iter_bobyqa = 3000
solver="bobyqa"
p = 10

#########
# Data:
#########
# Data for MSFT:
dat <- read.csv("data/realized_volatility_measures_MSFT_train.csv")
measures <- c("rv_med_5min","rvKernel_mth_5min","rbp_cov_5min","rv_med_10min","rvKernel_mth_10min")
realized_measures <- dat[,measures]

# Identify the rows having NA and having zero observations:
id <- which(realized_measures<=0,arr.ind = TRUE)
realized_measures[id] <-  unlist(lapply(id[,2], function(x) mean(realized_measures[,x]))) 

# Square root transformation of realized measures:
X <- sqrt(realized_measures)

# Checking Kendall's tau:
cor(X,method="kendall")

# Making the mean of each component as 1:
mean_X <- apply(X,2,mean)
a <- lapply(1:length(mean_X),function(i) X[,i]/mean_X[i])

X <- do.call(cbind,a)
colnames(X) <- measures

# Plot data matrix:
plot.ts(X)

#############################################################
# Step 2: Calculate MLE (By profile likelihood estimation):
#############################################################

# Dimension of data and order of lags:
d <- ncol(X)
p <- p # Max-lag
q <- 0 # For our project q=0

X_train <- X[1:(nrow(X)-1),]
X_test <- X[nrow(X),]

# Grid of beta values and initial A matrix:
all_initial_values <- init.val(X=X_train,maxlag=p,n_beta=20)
init_A_matrix <- all_initial_values$init_A_matrix
(beta_grid <- all_initial_values$beta_list)

#########################################
# Example to get an estimate of runtime:
#########################################
# Compute MLE using profile likelihood by gridding over beta values:
# system.time({mle_output <- mle(X=X_train,maxlag=p,init_A_matrix = init_A_matrix
#                                ,beta_pars = 4)})

beta_i_index <- expand.grid(i_index = 1:d,beta = beta_grid)
beta_grid_df <- data.frame(beta=beta_grid,model_id=1:length(beta_grid))
beta_i_index_df <- inner_join(beta_i_index,beta_grid_df)

library(doParallel)

# Set up parallel estimation:
no_cores <- detectCores()
cl <- makeCluster(no_cores - 1)
registerDoParallel(cl)
registerDoSNOW(cl)

# print out the progress for every iteration
progress <- function(n) cat(sprintf("task %d is complete/n", n))
opts <- list(progress=progress)

mle_output_list <- foreach(i=1:nrow(beta_i_index_df), .packages=c("ACDm","nloptr","Rsolnp"),
                           .export = ls(globalenv()),.verbose = TRUE,.options.snow = opts) %dopar% {
                             tryCatch({mle_ith_comp(X=X_train,i.index = beta_i_index_df$i_index[i],maxlag=p
                                                    ,init_A_matrix = init_A_matrix,beta_par=beta_i_index_df$beta[i]
                                                    ,max_iter_bobyqa = max_iter_bobyqa,solver=solver,eps=eps)}
                                      ,error=function(e){return(NA)})
                           }

stopCluster(cl)

# Binding all outputs:
all_mle_list <- list()
for(i in 1:length(mle_output_list)){
  mle_list <- mle_output_list[[i]]
  if(is.na(mle_list[1])==T){
    all_mle_list[[i]] <- NA
  }else{
    mle_df = cbind.data.frame(beta_i_index_df[i,],t(mle_list$A_i_estimate),convergence = mle_list$conv_i
                              ,n_iterations = mle_list$iter_i,loglik_i = mle_list$loglik_i)
    colnames(mle_df) <- c(colnames(beta_i_index_df),
                          paste("A",1:(p*d),sep="_"),"convergence","n_iterations"
                          ,"loglik_i")
    all_mle_list[[i]] <- mle_df
  }
}

# Row binding all mle results as data frame:
all_mle_df <- do.call(rbind,all_mle_list)

# Take only the complete cases:
all_mle_df <- all_mle_df[complete.cases(all_mle_df),]

# Filter out the models with incomplete information:
split_all_mle_df_bymodelid <- split.data.frame(all_mle_df,f = all_mle_df$model_id)
rm_id <- which(lapply(split_all_mle_df_bymodelid,nrow)<d)
if(length(rm_id)>0){
  split_all_mle_df_bymodelid <- split_all_mle_df_bymodelid[-rm_id]
}

# Adding full data log-lik in the data frame:
for(i in 1:length(split_all_mle_df_bymodelid)){
  split_all_mle_df_bymodelid[[i]]$loglik_full_data <- sum(split_all_mle_df_bymodelid[[i]]$loglik_i)
}

# Calculating log-likelihood by model_id:
mle_id <- which.max(unlist(lapply(1:length(split_all_mle_df_bymodelid),function(i) unique(split_all_mle_df_bymodelid[[i]]$loglik_full_data))))
A_matrix_estimated_mle <- split_all_mle_df_bymodelid[[mle_id]][,paste("A",1:(p*d),sep="_")]
beta_mle <- unique(split_all_mle_df_bymodelid[[mle_id]][,"beta"])

# A matrix list:
A_matrices <- A_matrix_estimated_mle
start_id <- seq(1,length(A_matrices),by=d)
end_id <- seq(d,length(A_matrices),by=d)
A_matrix_list_estimated <- list()

for(i in 1:length(start_id)){
  A_matrix_list_estimated[[i]] <- A_matrices[start_id[i]:end_id[i]]
}

# Writing omega in terms of the elements of matrix A and data
mu_vec <- apply(X=X_train,2,mean)
variance_vec <- diag(cov(X_train))
beta_par <- beta_mle
M <- rep(digamma(beta_par),d) - log(beta_par)
omega_mle <- crossprod(t((diag(x = 1,nrow = d,ncol = d) - Reduce(f = "+",x = A_matrix_list_estimated)))
                       , matrix((log(mu_vec) - variance_vec/(2*mu_vec^2)))) - M

##################################################
# Calculate lambda path (first get lambda_max):
##################################################
v <- 0.1 # If you want more sparsity, choose a small v
N <- nrow(X)
lambda_max <- abs(log_lik_logvMEM(X = X,A_matrix = A_matrix_estimated_mle
                                  ,beta_par = beta_mle))/(N*p*v)
epsilon <- 1e-2
K <- 20 # Number of lambda penalties
lambda_grid <- round(exp(seq(log(lambda_max), log(lambda_max*epsilon),
                            length.out = K)), digits = 10)

beta_lamba_combinations <- expand.grid(beta =beta_mle,lambda=lambda_grid)
beta_lamba_combinations$Model_id <- paste0("Model",1:nrow(beta_lamba_combinations))
beta_lamba_i_index_combinations <- expand.grid(i.index = 1:d, beta =beta_mle,lambda=lambda_grid)

#######################################################
# Parallel computing for getting parameter estimates:
#######################################################
library(doParallel)

# Combinations of penalty and HLag structure:
all_penalties <- c("group-lasso","group-mcp","group-scad","adaptive-group-lasso")
all_hlag <- c("componentwise","own-other","elementwise")

penalty_hlag_all <- expand.grid(penalty=all_penalties,Hlag_structure=all_hlag)

for(comb_id in 1:nrow(penalty_hlag_all)){
  Hlag_structure = as.character(penalty_hlag_all$Hlag_structure[comb_id])
  penalty = as.character(penalty_hlag_all$penalty[comb_id])
  threshold=threshold;maxlag <- p
  
  # Set up parallel estimation:
  no_cores <- detectCores()
  cl <- makeCluster(no_cores - 1)
  registerDoParallel(cl)
  registerDoSNOW(cl)
  
  # print out the progress for every iteration
  progress <- function(n) cat(sprintf("task %d is complete/n", n))
  opts <- list(progress=progress)
  
  output_list <- foreach(beta_lambda_index=1:nrow(beta_lamba_i_index_combinations), .packages=c("ACDm","nloptr","Rsolnp"),
                         .export = ls(globalenv()),.verbose = TRUE,.options.snow = opts) %dopar% {
                           tryCatch({penalized_logvMEM_fix_beta_fix_lambda(X_train,X_test,i.index=beta_lamba_i_index_combinations$i.index[beta_lambda_index]
                                                                           ,maxlag=maxlag,max_iter_bobyqa=max_iter_bobyqa,eps = eps,penalty=penalty
                                                                           ,Hlag_structure=Hlag_structure,verbose=TRUE
                                                                           ,lambda_penalty = beta_lamba_i_index_combinations$lambda[beta_lambda_index]
                                                                           ,init_A_i_list=lapply(1:d,function(i.index) all_initial_values$init_A_matrix[i.index,])
                                                                           ,threshold=threshold,mle_A_matrix=A_matrix_estimated_mle
                                                                           ,beta_par = beta_lamba_i_index_combinations$beta[beta_lambda_index]
                                                                           ,h_step=1,solver=solver)},error=function(e){return(NA)
                                                                           })
                         }
  # Stop clusters:
  stopCluster(cl)
  
  # Binding all outputs:
  combined_output_df <- do.call(rbind,output_list)
  library(dplyr)
  output_table <- inner_join(x=combined_output_df,y = beta_lamba_combinations)
  
  # Final output table:
  output_table <- output_table[complete.cases(output_table),]
  
  # Split the output table by model id:
  split_output_bymodel <- split.data.frame(output_table,f = output_table$Model_id)
  
  # Check which model has complete information:
  rm_id <- which(lapply(1:length(split_output_bymodel),function(x) nrow(split_output_bymodel[[x]]))<d)
  
  if(length(rm_id) >0){
    split_output_bymodel <- ifelse(length(rm_id) > 0,split_output_bymodel[-rm_id],split_output_bymodel)
  }
  
  ###############################
  # Compute full data metrics:
  ###############################
  model_result_list <- list()
  for(i in 1:length(split_output_bymodel)){
    model_data <- split_output_bymodel[[i]]  
    A_matrix_estimated <- model_data[,paste("A",1:(d*p),sep="_")]
    beta_par <- unique(model_data[,"beta"])
    model_data$loglik_full_data <- log_lik_logvMEM(X = X,A_matrix = A_matrix_estimated,beta_par = beta_par)
    
    # AIC, BIC, MAD (in-sample), df and penalized log-lik:
    model_data$df_full_data <- nrow(which(abs(A_matrix_estimated)>threshold,arr.ind = T))
    model_data$aic_full_data <- 2*model_data$df_full_data - 2*model_data$loglik_full_data
    model_data$bic_full_data <- model_data$df_full_data*log(N) - 2*model_data$loglik_full_data 
    model_data$MAD_in_sample <- mean(model_data$MAD_i)
    model_data$count_non_zero_coef_full_data=sum(model_data$count_non_zero_coef_i)
    model_data$pen_loglik_full_data <- sum(model_data$pen_loglik_i)
    omega_vec <- model_data$omega_i
    
    # h-step forecast metrics:
    forecast_metrics_logvMEM_obj <- forecast_logvMEM_full_data(X_train=X_train,X_test=X_test,h_step=1
                                                                     ,omega_estimated=omega_vec,A_matrix_estimated=A_matrix_estimated,maxlag=maxlag)
    forecast_metrics_logvMEM_i <- forecast_metrics_logvMEM_obj$forecast_df
    model_data = cbind.data.frame(model_data,forecast_metrics_logvMEM_i[,c("MSFE_i","MAD_i","MAPE_i")])
    
    forecast_metrics_logvMEM_full_data <- forecast_metrics_logvMEM_obj$model_metrics
    model_data$MSFE_forecast <- forecast_metrics_logvMEM_full_data["MSFE"]
    model_data$MAD_forecast <- forecast_metrics_logvMEM_full_data["MAD"]
    model_data$MAPE_forecast <- forecast_metrics_logvMEM_full_data["MAPE"]
    model_result_list[[i]] <- model_data
  }
  
  result_df <- do.call(rbind,model_result_list)
  
  # Setting output directory:
  output_dir = paste0(project_dir,"real_data_analysis/output_wd_sqrt_transformation/")
  setwd(output_dir)
  
  # Write the output table to a csv file:
  file_name_output <- paste0(paste("output_table",Hlag_structure,penalty,sep="-"),".csv")
  write.csv(result_df,file_name_output,row.names = F)
  
  # Split result_df by lambda:
  split_result_df_by_model_comp <- split.data.frame(result_df,f = result_df$component)
  
  # Choose the lambda for each component with componentwise min BIC:
  final_estimate_df_bic=do.call(rbind,lapply(1:length(split_result_df_by_model_comp)
                                             ,function(x) split_result_df_by_model_comp[[x]][which.min(split_result_df_by_model_comp[[x]]$bic_i),]))
  
  # Write the final estimate to a csv file:
  file_name_final_estimate <- paste0(paste("final_estimate_bic",Hlag_structure,penalty,sep="-"),".csv")
  write.csv(final_estimate_df_bic,file_name_final_estimate,row.names = F)
  
  # Choose the lambda for each component with componentwise min MSFE:
  final_estimate_df_msfe=do.call(rbind,lapply(1:length(split_result_df_by_model_comp)
                                             ,function(x) split_result_df_by_model_comp[[x]][which.min(split_result_df_by_model_comp[[x]]$MSFE_i),]))
  
  # Write the final estimate to a csv file:
  file_name_final_estimate <- paste0(paste("final_estimate_MSFE",Hlag_structure,penalty,sep="-"),".csv")
  write.csv(final_estimate_df_msfe,file_name_final_estimate,row.names = F)
  
  end_time <- Sys.time()
  
  paste("Start time is ",start_time)
  paste("End time is ",end_time)
  
}










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
source("/Users/madhusreechowdhury/Desktop/Project_2_Revised 2/codes/sparse_logvMEM_estimation.R")
source("/Users/madhusreechowdhury/Desktop/Project_2_Revised 2/codes/sparsity_plot.R")

# Global parameters that control the optimization:
eps = 1e-3;threshold = 1e-3
max_iter_bobyqa = 3000
solver="bobyqa"
Hlag_structure = "componentwise"

# Setting working directory:
setwd("/Users/madhusreechowdhury/Desktop/Project_2_Revised 2/simulation_study_v3/simulation_study_output/componentwise/")

# Setting random seed for reproducibility:
set.seed(12345)

#####################################
# Step 1: Generating simulated data:
#####################################
d <- 5
p <- 10
q <- 0
r <- max(p,q)

# Max-lag matrices:
L <- matrix(c(rep(3,d),rep(4,d),rep(4,d),rep(5,d),rep(5,d)),nrow=d,ncol=d,byrow = TRUE) # HLag- componentwise

# List of A matrices:
A_matrix_list <- autoregressive_matrix(d = d,p = p,L = L,max_eigen_abs = 0.96)$A_matrix_list

# A matrix:
A_matrix <- do.call(cbind,A_matrix_list)

# Writing A matrix:
write.csv(A_matrix,paste0(paste("A_matrix",Hlag_structure,sep="_"),".csv"),row.names = F)

# Data: (Generating 1500 train data points and 1 for testing) 
X <- sim_logvMEM(N = 1501,d = d,omega=matrix(rep(0.1,d),ncol=1),A=A_matrix_list
                 ,B=NULL,beta_par=4,burn_in = 500)

# Writing data matrix:
write.csv(X,paste0(paste("data_matrix",Hlag_structure,sep="_"),".csv"),row.names = F)

# Summary of X:
apply(X,2,mean)

# Covariance matrix:
var(X)

# Correlation matrix:
cor(X,method = "kendall")

#############################################################
# Step 2: Calculate MLE (By profile likelihood estimation):
#############################################################
p <- p # Max-lag
q <- 0 # For our project q=0

X_train <- X[1:1500,]
X_test <- X[1501,]

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
progress <- function(n) cat(sprintf("task %d is complete\n", n))
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

# Compare the actual and the mle estimates:
compare_df <- data.frame(actual = unlist(lapply(1:nrow(A_matrix),function(i) A_matrix[i,]))
                         , pred = unlist(lapply(1:nrow(A_matrix_estimated_mle),function(i) A_matrix_estimated_mle[i,])))

# Writing omega in terms of the elements of matrix A and data
mu_vec <- apply(X=X_train,2,mean)
variance_vec <- diag(cov(X_train))
beta_par <- beta_mle
M <- rep(digamma(beta_par),d) - log(beta_par)
omega_mle <- crossprod(t((diag(x = 1,nrow = d,ncol = d) - Reduce(f = "+",x = A_matrix_list)))
                       , matrix((log(mu_vec) - variance_vec/(2*mu_vec^2)))) - M

##################################################
# Calculate lambda path (first get lambda_max):
##################################################
v <- 0.1 # If you want more sparsity, choose a small v
N <- nrow(X)
lambda_max <- abs(log_lik_logvMEM(X = X,A_matrix = A_matrix_estimated_mle
                                  ,beta_par = beta_mle))/(N*p*v)
epsilon <- 1e-3
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

all_penalties = c("group-lasso","group-scad","group-mcp")

for(penalty in all_penalties){
  #Hlag_structure="componentwise"
  threshold=threshold;maxlag <- p
  
  # Set up parallel estimation:
  no_cores <- detectCores()
  cl <- makeCluster(no_cores - 1)
  registerDoParallel(cl)
  registerDoSNOW(cl)
  
  # print out the progress for every iteration
  progress <- function(n) cat(sprintf("task %d is complete\n", n))
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
  
  # Setting working directory:
  setwd("/Users/madhusreechowdhury/Desktop/Project_2_Revised 2/simulation_study_v3/simulation_study_output/componentwise/")
  
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
  
  # Sparsity plot of actual A matrix:
  A_matrix_sparse = A_matrix
  A_matrix_sparse[which(abs(A_matrix_sparse)<=threshold)]=0
  A_matrix_sparse[which(abs(A_matrix_sparse)>threshold)]=1
  
  filename1 = paste(paste(Hlag_structure,"actual_plot",sep="_"),"pdf",sep=".")
  pdf(filename1)
  SparsityPlot(B = A_matrix_sparse,p = p,k = d,s = 2,m = 2,title = filename1)
  dev.off()
  
  # Sparsity plot of the estimated final matrix based on BIC:
  A_penalized_estimate = final_estimate_df_bic[,paste("A",1:(p*d),sep="_")]
  A_penalized_estimate[which(abs(A_penalized_estimate)<=threshold,arr.ind = T)]=0
  A_penalized_estimate[which(abs(A_penalized_estimate)>threshold,arr.ind = T)]=1
  
  filename2 = paste(paste(penalty,Hlag_structure,"estimated_plot_bic",sep="_"),"pdf",sep=".")
  pdf(filename2)
  SparsityPlot(B = A_penalized_estimate,p = p,k = d,s = 2,m = 2,title = filename2)
  dev.off()
  
  # Sparsity plot of the estimated final matrix based on MSFE:
  A_penalized_estimate = final_estimate_df_msfe[,paste("A",1:(p*d),sep="_")]
  A_penalized_estimate[which(abs(A_penalized_estimate)<=threshold,arr.ind = T)]=0
  A_penalized_estimate[which(abs(A_penalized_estimate)>threshold,arr.ind = T)]=1
  
  filename3 = paste(paste(penalty,Hlag_structure,"estimated_plot_MSFE",sep="_"),"pdf",sep=".")
  pdf(filename3)
  SparsityPlot(B = A_penalized_estimate,p = p,k = d,s = 2,m = 2,title = filename3)
  dev.off()
  
}








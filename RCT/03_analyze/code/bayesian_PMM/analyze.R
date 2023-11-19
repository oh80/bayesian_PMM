main <- function(){
  #read data
  missing_type <- "NMAR"
  missing_rate <- 0.3
  
  file_name <- paste0(missing_type,"_", missing_rate, ".", "obj")
  path <- here::here("RCT","02_build","data", missing_type ,file_name)
  data <- readRDS(path)
  
  #read functions
  code_path <- here::here("RCT","03_analyze","code","bayesian_PMM","gibbs_sampler.R")
  source(code_path)
  code_path <- here::here("RCT","03_analyze","code","bayesian_PMM","get_prior_params.R")
  source(code_path)
  
  #set sample_size and burn_in
  sample_size <- 2000
  burn_in <- 400

  #complete data giggs sampler
  R3_params_sample <- R3_gibbs_sampler(data, sample_size) |> cut_burn_in(burn_in)
  R3_estimated_values <- get_estimated_values(R3_params_sample)
  
  #R=2 data gibbs sampler
  R2_prior_parameter <- get_R2_prior_parameter(R3_params_sample)
  R2_sample <- R2_gibbs_sampler(data, sample_size, R2_prior_parameter) |> cut_burn_in(burn_in)
  R2_estimated_values <- get_estimated_values(R2_sample)
  
  #R=1 data gibbs sampler
  R1_prior_parameter <- get_R1_prior_parameter(R2_sample, R3_params_sample)
  R1_sample <- R1_gibbs_sampler(data, sample_size, R1_prior_parameter) |> cut_burn_in(burn_in)
  R1_estimated_values <- get_estimated_values(R1_sample)
  
  #R=0 data gibbs sampler
  R0_prior_parameter <- get_R0_prior_parameter(R1_sample, R2_sample, R3_params_sample)
  R0_sample <- R0_gibbs_sampler(data, sample_size, R0_prior_parameter) |> cut_burn_in(burn_in)
  R0_estimated_values <- get_estimated_values(R0_sample)
  
  #combine results
  standard_error <- get_se(R2_prior_parameter, R1_prior_parameter, R0_prior_parameter,
                           R3_params_sample, R2_sample, R1_sample, R0_sample, data)
  results <- combine_results(R3_estimated_values, R2_estimated_values,
                             R1_estimated_values, R0_estimated_values,
                             data, standard_error)
  
  #save estimated values
  file_name <- paste0(missing_type ,missing_rate,".obj")
  #save(estimate_result, file_name)
  return(results)
}


get_estimated_values <- function(sample){
  N <- length(sample$beta1)
  beta1 <- sample$beta1[[1]]/N
  for (i in 2:N){
    add <- sample$beta1[[i]]/N
    beta1 <- beta1+add
  }
  
  N <- length(sample$beta2)
  beta2 <- sample$beta2[[1]]/N
  for (i in 2:N){
    add <- sample$beta2[[i]]/N
    beta2 <- beta2+add
  }
  
  N <- length(sample$beta3)
  beta3 <- sample$beta3[[1]]/N
  for (i in 2:N){
    add <- sample$beta3[[i]]/N
    beta3 <- beta3+add
  }
  
  output <- list()
  output$beta1 <- beta1
  output$beta2 <- beta2
  output$beta3 <- beta3
  return(output)
}


cut_burn_in <- function(sample, burn_in){
  sample_size <- length(sample$beta1)
  variable_name <- names(sample)
  output <- list()
  
  for (k in 1:length(sample)) {
    output[[k]] <- sample[[k]][burn_in:sample_size]
  }
  
  names(output) <- variable_name
  return(output)
}


get_se <- function(R2_prior_parameter, R1_prior_parameter, R0_prior_parameter,
                   R3_sample, R2_sample,R1_sample,  R0_sample, data){
  n <- data |> dplyr::group_by(R)|> dplyr::count(name = "R")
  n0 <- n[1,1] |> as.numeric()
  n1 <- n[2,1] |> as.numeric()
  n2 <- n[3,1] |> as.numeric()
  n3 <- n[4,1] |> as.numeric()
  N <- length(data$R)
  
  #Sigma_3
  R3_sigma_3 <- R2_prior_parameter[[6]]
  R2_sigma_3 <- R1_prior_parameter[[6]]
  R1_sigma_3 <- R0_prior_parameter[[6]]
  
  R0_mu_3 <- R0_sample$beta3[[1]] / length(R0_sample$beta3)
  for (i in 2:length(R0_sample$beta3)){
    add <- R0_sample$beta3[[i]] / length(R0_sample$beta3)
    R0_mu_3 <- R0_mu_3 + add
  }
  
  R0_sigma_3 <- t(R0_sample$beta3[[1]] - R0_mu_3)  %*% (R0_sample$beta3[[1]] - R0_mu_3) / length(R0_sample$beta3)
  for (i in 2:length(R0_sample$beta3)){
    add <- t(R0_sample$beta3[[i]] - R0_mu_3)  %*% (R0_sample$beta3[[i]] - R0_mu_3) / length(R0_sample$beta3)
    R0_sigma_3<- R0_sigma_3 + add
  }
  
  sigma_3 <- R3_sigma_3 * n3/N + R2_sigma_3  * n2/N + R1_sigma_3  * n1/N + R0_sigma_3 * n0/N
  
  #Sigma_2
  R2_sigma_2 <- R1_prior_parameter[[4]]
  R1_sigma_2 <- R0_prior_parameter[[4]]
  
  R3_mu_2 <- R3_sample$beta2[[1]] / length(R3_sample$beta2)
  for (i in 2:length(R3_sample$beta2)){
    add <- R3_sample$beta2[[i]] / length(R3_sample$beta2)
    R3_mu_2 <- R3_mu_2 + add
  }
  
  R3_sigma_2 <- t(R3_sample$beta2[[1]] - R3_mu_2)  %*% (R3_sample$beta2[[1]] - R3_mu_2) / length(R3_sample$beta2)
  for (i in 2:length(R0_sample$beta2)){
    add <- t(R3_sample$beta2[[i]] - R3_mu_2)  %*% (R3_sample$beta2[[i]] - R3_mu_2) / length(R3_sample$beta2)
    R3_sigma_2<- R3_sigma_2 + add
  }
  
  R0_mu_2 <- R0_sample$beta2[[1]] / length(R0_sample$beta2)
  for (i in 2:length(R0_sample$beta2)){
    add <- R0_sample$beta2[[i]] / length(R0_sample$beta2)
    R0_mu_2 <- R0_mu_2 + add
  }
  
  R0_sigma_2 <- t(R0_sample$beta2[[1]] - R0_mu_2)  %*% (R0_sample$beta2[[1]] - R0_mu_2) / length(R0_sample$beta2)
  for (i in 2:length(R0_sample$beta2)){
    add <- t(R0_sample$beta2[[i]] - R0_mu_2)  %*% (R0_sample$beta2[[i]] - R0_mu_2) / length(R0_sample$beta2)
    R0_sigma_2<- R0_sigma_2 + add
  }
  
  sigma_2 <- R3_sigma_2 * n3/N + R2_sigma_2  * n2/N + R1_sigma_2  * n1/N + R0_sigma_2 * n0/N
  
  #Sigma_1
  R1_sigma_1 <- R0_prior_parameter[[2]]
  R3_mu_1 <- R3_sample$beta1[[1]] / length(R3_sample$beta1)
  for (i in 2:length(R3_sample$beta1)){
    add <- R3_sample$beta1[[i]] / length(R3_sample$beta1)
    R3_mu_1 <- R3_mu_1 + add
  }
  
  R3_sigma_1 <- t(R3_sample$beta1[[1]] - R3_mu_1)  %*% (R3_sample$beta1[[1]] - R3_mu_1) / length(R3_sample$beta1)
  for (i in 2:length(R3_sample$beta2)){
    add <- t(R3_sample$beta1[[i]] - R3_mu_1)  %*% (R3_sample$beta1[[i]] - R3_mu_1) / length(R3_sample$beta1)
    R3_sigma_1<- R3_sigma_1 + add
  }
  
  R2_mu_1 <- R2_sample$beta1[[1]] / length(R2_sample$beta1)
  for (i in 2:length(R2_sample$beta1)){
    add <- R2_sample$beta1[[i]] / length(R2_sample$beta1)
    R2_mu_1 <- R2_mu_1 + add
  }
  
  R2_sigma_1 <- t(R2_sample$beta1[[1]] - R2_mu_1)  %*% (R2_sample$beta1[[1]] - R2_mu_1) / length(R2_sample$beta1)
  for (i in 2:length(R2_sample$beta2)){
    add <- t(R2_sample$beta1[[i]] - R2_mu_1)  %*% (R2_sample$beta1[[i]] - R2_mu_1) / length(R2_sample$beta1)
    R2_sigma_1<- R2_sigma_1 + add
  }
  
  R0_mu_1 <- R0_sample$beta1[[1]] / length(R0_sample$beta1)
  for (i in 2:length(R0_sample$beta1)){
    add <- R0_sample$beta1[[i]] / length(R0_sample$beta1)
    R0_mu_1 <- R0_mu_1 + add
  }
  
  R0_sigma_1 <- t(R0_sample$beta1[[1]] - R0_mu_1)  %*% (R0_sample$beta1[[1]] - R0_mu_1) / length(R0_sample$beta1)
  for (i in 2:length(R0_sample$beta2)){
    add <- t(R0_sample$beta1[[i]] - R0_mu_1)  %*% (R0_sample$beta1[[i]] - R0_mu_1) / length(R0_sample$beta1)
    R0_sigma_1<- R0_sigma_1 + add
  }
  
  sigma_1 <- R3_sigma_1 * n3/N + R2_sigma_1  * n2/N + R1_sigma_1  * n1/N + R0_sigma_1 * n0/N
  
  output <- list()
  output$sigma_3 <- sigma_3
  output$sigma_2 <- sigma_2
  output$sigma_1 <- sigma_1
  return(output)
}


combine_results <- function(R3_estimated_values, R2_estimated_values,
                            R1_estimated_values, R0_estimated_values,
                            data, standard_error){
  n <- data |> dplyr::group_by(R)|> dplyr::count(name = "R")
  n0 <- n[1,1] |> as.numeric()
  n1 <- n[2,1] |> as.numeric()
  n2 <- n[3,1] |> as.numeric()
  n3 <- n[4,1] |> as.numeric()
  N <- length(data$R)
  
  beta1 <- R3_estimated_values$beta1 * n3/N + R2_estimated_values$beta1 * n2/N +
    R1_estimated_values$beta1 * n1/N + R0_estimated_values$beta1 * n0/N
  
  beta2 <- R3_estimated_values$beta2 * n3/N + R2_estimated_values$beta2 * n2/N +
    R1_estimated_values$beta2 * n1/N + R0_estimated_values$beta2 * n0/N
  
  beta3 <- R3_estimated_values$beta3 * n3/N + R2_estimated_values$beta3 * n2/N +
    R1_estimated_values$beta3 * n1/N + R0_estimated_values$beta3 * n0/N
  
  output <- list()
  output$beta1 <- beta1
  output$beta2 <- beta2
  output$beta3 <- beta3
  output$standard_error <- standard_error
  return(output)
}


save <- function(results, file_name){
  path <- here::here("RCT", "03_analyze","output", "bayesian_PMM", file_name)
  saveRDS(results, path)
}



results <- main()

results$beta1
results$beta2
results$beta3

results$standard_error$sigma_3
results$standard_error$sigma_2
results$standard_error$sigma_1
main <- function(){
  #read functions
  code_path <- here::here("Liner_Regression","03_analyze","code","Bayesian_PMM","compute_complete_data_params.R")
  source(code_path)
  code_path <- here::here("Liner_Regression","03_analyze","code","Bayesian_PMM","gibbs_sampler.R")
  source(code_path)
  
  # read data 
  missing_rate <- 0.5
  missing_type <- "NMAR"
  file_name <- paste0(missing_type , "_", missing_rate, ".obj")
  path <- here::here("Liner_regression", "02_build","data", file_name)
  data <- readRDS(path)
  
  # set hyper parameter
  g <- 1
  sigma <- get_sigma(data)
  sample_size <- 1000
  
  beta_0 <- get_OLS_extimater(data)
  sigma_0 <- get_Sigma0(data, g, sigma)
  
  # generate sample
  sample <- generate_sample(data, sample_size, beta_0, sigma_0, sigma)
  
  # save sample
  save(sample, file_name)
}


generate_sample <- function(data, sample_size, beta_0, Sigma_0, sigma){
  output <- list()
  for (i in 1:10) {
    sample <- gibbs_sampler(data, j = i, sample_size, beta_0, Sigma_0, sigma)
    output[[i]] <- sample
  }
  return(output)
}


save <- function(sample, file_name){
  path <- here::here("Liner_regression", "03_analyze","output", "bayesian_PMM", file_name)
  saveRDS(sample, path)
}


main()


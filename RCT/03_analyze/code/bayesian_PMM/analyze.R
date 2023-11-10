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

  #get complete data giggs sampler
  R3_params_sample <- R3_gibbs_sampler(data, 100)
  
  return(R3_params_sample)
  
}

sample <- main()
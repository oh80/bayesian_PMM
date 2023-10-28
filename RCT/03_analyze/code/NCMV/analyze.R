main <- function(){
  #read data
  missing_type <- "MCAR"
  missing_rate <- 0.1
  
  file_name <- paste0(missing_type,"_", missing_rate, ".", "obj")
  path <- here::here("RCT","02_build","data", missing_type ,file_name)
  data <- readRDS(path)
  
  #read functions
  code_path <- here::here("RCT","03_analyze","code","NCMV","imputation.R")
  source(code_path)
  
  #param set
  D <- 10
  delta <- 0
  
  #imputation
  psude_complete_data <- multiple_imputation(data, D, delta)
  
  #estimation
  estimated_values <- esimate_each_data(psude_complete_data)
  
  #combine results
  estimate_result <- combine(estimated_values)

  #save estimated values
  file_name <- paste0(missing_type ,missing_rate,"delta",delta,".obj")
  save(estimate_result, file_name)
}


multiple_imputation <- function(data, D, delta){
  output <- list()
  weight1 <- c(0,0,1)
  weight2 <- c(0,1)
  for (d in 1:D) {
    set.seed(d)
    output[[d]] <- single_imputation(data, delta)
  }
  return(output)
}


esimate_each_data <- function(psude_complete_data){
  D <- length(psude_complete_data)
  output <- list()
  for (d in 1:D) {
    data <- psude_complete_data[[d]]
    estimated_values <- list()
    estimated_values[[1]] <- estimatr::lm_robust(formula = Y1 ~ X1+X2+X3+X4+X5+Treatment,
                                                 data = data,se_type = "stata")
    estimated_values[[2]] <- estimatr::lm_robust(formula = Y2 ~ X1+X2+X3+X4+X5+Treatment,
                                                 data = data,se_type = "stata")
    estimated_values[[3]] <- estimatr::lm_robust(formula = Y3 ~ X1+X2+X3+X4+X5+Treatment,
                                                 data = data,se_type = "stata")
    output[[d]] <- estimated_values
    
  }
  return(output)
}


combine <- function(estimated_values){
  D = length(estimated_values)
  
  combined_values1 <- rep(0, 7)
  for (d in 1:D) {
    combined_values1 <- combined_values1 + estimated_values[[d]][[1]]$coefficients/D
  }
  
  combined_values2 <- rep(0, 7)
  for (d in 1:D) {
    combined_values2 <- combined_values2 + estimated_values[[d]][[2]]$coefficients/D
  }
  
  combined_values3 <- rep(0, 7)
  for (d in 1:D) {
    combined_values3 <- combined_values3 + estimated_values[[d]][[3]]$coefficients/D
  }
  
  var1 <- rep(0, 7)
  for (d in 1:D) {
    var1 <- var1 + 1/D * (estimated_values[[d]][[1]]$std.error)^2 + 
      (1+1/D)*(1/(D-1))*(combined_values1 - estimated_values[[d]][[1]]$coefficients)^2
  }
  
  var2 <- rep(0, 7)
  for (d in 1:D) {
    var2 <- var2 + 1/D * (estimated_values[[d]][[2]]$std.error)^2 + 
      (1+1/D)*(1/(D-1))*(combined_values2 - estimated_values[[d]][[2]]$coefficients)^2
  }
  
  var3 <- rep(0, 7)
  for (d in 1:D) {
    var3 <- var3 + 1/D * (estimated_values[[d]][[3]]$std.error)^2 + 
      (1+1/D)*(1/(D-1))*(combined_values3 - estimated_values[[d]][[3]]$coefficients)^2
  }
  
  output <- list()
  output$estimated_values1 <- combined_values1
  output$estimated_values2 <- combined_values2
  output$estimated_values3 <- combined_values3
  output$var1 <- var1
  output$var2 <- var2
  output$var3 <- var3
  
  return(output)
}


save <- function(results, file_name){
  path <- here::here("RCT", "03_analyze","output", "NCMV", file_name)
  saveRDS(results, path)
}


main()



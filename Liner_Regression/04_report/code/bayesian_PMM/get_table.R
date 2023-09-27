main <- function(){
  code_path <- here::here("Liner_Regression","04_report","code","bayesian_PMM","get_estimated_values.R")
  source(code_path)
  
  #set params
  missing_rate <- 0.1
  
  #read results
  results <- read_results(missing_rate)
  sd <- read_sd(missing_rate)
  
  #get table
  table <- get_data_frame(results, sd) |> get_table()
  
  #save table as text
  save(table, missing_rate)
}


read_results <- function(missing_rate){
  estimated_values1 <- get_estimated_values(missing_rate, missing_type = "MCAR")
  estimated_values2 <- get_estimated_values(missing_rate, missing_type = "MAR")
  estimated_values3 <- get_estimated_values(missing_rate, missing_type = "NMAR")
  
  output <- list()
  output[[1]] <- estimated_values1
  output[[2]] <- estimated_values2
  output[[3]] <- estimated_values3
  return(output)
}

read_sd <- function(missing_rate){
  estimated_values1 <- get_standard_error(missing_rate, missing_type = "MCAR")
  estimated_values2 <- get_standard_error(missing_rate, missing_type = "MAR")
  estimated_values3 <- get_standard_error(missing_rate, missing_type = "NMAR")
  
  output <- list()
  output[[1]] <- estimated_values1
  output[[2]] <- estimated_values2
  output[[3]] <- estimated_values3
  return(output)
}


get_data_frame <- function(results, sd){
  effect_col <- paste0("beta",seq(1, 10, by=1))
  MCAR <- as.vector(results[[1]])
  MAR <- as.vector(results[[2]])
  NMAR <- as.vector(results[[3]])
  
  coef_df <- data.frame(effect = effect_col,
                        MCAR = sprintf("%.1f",MCAR),
                        MAR = sprintf("%.1f",MAR),
                        NMAR = sprintf("%.1f",NMAR))
  
  MCAR_sd = as.vector(t(sd[[1]]))
  MAR_sd = as.vector(t(sd[[2]]))
  NMAR_sd = as.vector(t(sd[[3]]))
  
  sd_df <- data.frame(effect = effect_col,
                      MCAR = format(MCAR_sd, digits = 2),
                      MAR = format(MAR_sd, digits = 2),
                      NMAR = format(NMAR_sd, digits = 2))
  
  output <- data.frame(effect = coef_df$effect,
                       MCAR = stringr::str_c(coef_df[[2]], " (",sd_df[[2]], ")"),
                       MAR = stringr::str_c(coef_df[[3]], " (",sd_df[[3]], ")"),
                       NMAR = stringr::str_c(coef_df[[4]], " (",sd_df[[4]], ")"))
  
  return(output)
}


get_table <- function(data_frame){
  output <- data_frame |> kableExtra::kable(booktabs = TRUE, format = "latex")
  
  return(output)
}



save <- function(table, missing_rate){
  file_name <- paste0(missing_rate,"_table.txt")
  path <- here::here("Liner_Regression","04_report","output","bayesian_PMM", file_name)
  writeLines(table, path)
}


main()




main <- function(){
  #set params
  missing_rate <- 0.1
  
  #read results
  results <- read_results(missing_rate)
  
  #get table
  table <- results |> get_data_frame() |> 
    get_table()
  
  #save table as text
  save(table, missing_rate)
}


read_results <- function(missing_rate){
  file_name1 <- paste0("MCAR_",missing_rate,".obj")
  file_name2 <- paste0("MAR_",missing_rate,".obj")
  file_name3 <- paste0("NMAR_",missing_rate,".obj")
  
  result1 <- readRDS(here::here("Liner_Regression","03_analyze","output","CCMV",file_name1))
  result2 <- readRDS(here::here("Liner_Regression","03_analyze","output","CCMV",file_name2))
  result3 <- readRDS(here::here("Liner_Regression","03_analyze","output","CCMV",file_name3))
  
  output <- list()
  output[[1]] <- result1
  output[[2]] <- result2
  output[[3]] <- result3
  return(output)
}


get_data_frame <- function(results){
  effect_col <- paste0("beta",seq(1, 10, by=1))
  MCAR <- as.vector(results[[1]]$estimated_values[2:11])
  MAR <- as.vector(results[[2]]$estimated_values[2:11])
  NMAR <- as.vector(results[[3]]$estimated_values[2:11])
  
  coef_df <- data.frame(effect = effect_col,
                        MCAR = format(MCAR, digits = 1),
                        MAR = format(MAR, digits = 1),
                        NMAR = format(NMAR, digits = 1))
  
  MCAR_sd = sqrt(results[[1]]$var[2:11])
  MAR_sd = sqrt(results[[2]]$var[2:11])
  NMAR_sd = sqrt(results[[2]]$var[2:11])
  
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
  path <- here::here("Liner_Regression","04_report","output","CCMV", file_name)
  writeLines(table, path)
}


main()



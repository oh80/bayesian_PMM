main <- function(){
  library(ggplot2)
  
  #read missing 
  missing_rate <- 0.3
  results <- read_results(missing_rate)
  
  
  #get table
  results <- read_results(missing_rate) 
  table <- results |> get_data_frame() |> get_table()
  
  #save table as text
  save(table, missing_rate)
}


read_results <- function(missing_rate){
  path1 <- here::here("RCT", "03_analyze", "output", "bayesian_PMM",
                      paste0("MCAR", missing_rate,".obj" ))
  path2 <- here::here("RCT", "03_analyze", "output", "bayesian_PMM",
                      paste0("MAR", missing_rate,".obj" ))
  path3 <- here::here("RCT", "03_analyze", "output", "bayesian_PMM",
                      paste0("NMAR", missing_rate,".obj" ))
  
  results1 <- readRDS(path1)
  results2 <- readRDS(path2)
  results3 <- readRDS(path3)

  output <- list(results1, results2, results3)
  return(output)
}

get_data_frame <- function(results){
  time <- c(paste0(seq(1, 3, by=1)))
  MCAR <- c(results[[1]]$beta1[2], results[[1]]$beta2[2], results[[1]]$beta3[2])
  MAR <- c(results[[2]]$beta1[2], results[[2]]$beta2[2], results[[2]]$beta3[2])
  NMAR <- c(results[[3]]$beta1[2], results[[3]]$beta2[2], results[[3]]$beta3[2])
  
  MCAR_sd <- c(results[[1]]$standard_error$sigma_1[2,2], results[[1]]$standard_error$sigma_2[2,2],
               results[[1]]$standard_error$sigma_3[2,2])
  MAR_sd <- c(results[[2]]$standard_error$sigma_1[2,2], results[[2]]$standard_error$sigma_2[2,2],
               results[[2]]$standard_error$sigma_3[2,2])
  NMAR_sd <- c(results[[3]]$standard_error$sigma_1[2,2], results[[3]]$standard_error$sigma_2[2,2],
               results[[3]]$standard_error$sigma_3[2,2])
  
  effect_df <- data.frame(
                          time = time,
                          MCAR = format(MCAR, digits = 2),
                          MAR = format(MAR, digits = 2),
                          NMAR = format(NMAR, digits = 2))
  sd_df <- data.frame(
                      time = time,
                      MCAR_sd = format(MCAR_sd, digits = 2),
                      MAR_sd = format(MAR_sd, digits = 2),
                      NMAR_sd = format(NMAR_sd, digits = 2))
  
  output <- data.frame(
                       time = time,
                       MCAR = stringr::str_c(effect_df[[2]], " (",sd_df[[2]], ")"),
                       MAR = stringr::str_c(effect_df[[3]], " (",sd_df[[3]], ")"),
                       NMAR = stringr::str_c(effect_df[[4]], " (",sd_df[[4]], ")"))
  return(output)
}


get_table <- function(data_frame){
  output <- data_frame |> kableExtra::kable(booktabs = TRUE, format = "latex")
  
  return(output)
}


save <- function(table, missing_rate){
  file_name <- paste0(missing_rate,"_table.txt")
  path <- here::here("RCT","04_report","output","bayesian_PMM", file_name)
  writeLines(table, path)
}


main()







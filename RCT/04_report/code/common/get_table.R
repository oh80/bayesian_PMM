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
  path1 <- here::here("RCT", "03_analyze", "output", "NCMV",
                      paste0("MCAR", missing_rate,"delta", 0,".obj" ))
  path2 <- here::here("RCT", "03_analyze", "output", "NCMV",
                      paste0("MAR", missing_rate,"delta", 0, ".obj" ))
  path3 <- here::here("RCT", "03_analyze", "output", "NCMV",
                      paste0("NMAR", missing_rate,"delta", 0 ,".obj" ))
  
  path4 <- here::here("RCT", "03_analyze", "output", "ACMV",
                      paste0("MCAR", missing_rate,"delta", 0 ,".obj" ))
  path5 <- here::here("RCT", "03_analyze", "output", "ACMV",
                      paste0("MAR", missing_rate,"delta", 0 ,".obj" ))
  path6 <- here::here("RCT", "03_analyze", "output", "ACMV",
                      paste0("NMAR", missing_rate,"delta",0 ,".obj" ))
  
  results1 <- readRDS(path1)
  results2 <- readRDS(path2)
  results3 <- readRDS(path3)
  
  results4 <- readRDS(path4)
  results5 <- readRDS(path5)
  results6 <- readRDS(path6)
  
  output <- list(results1, results2, results3, results4, results5, results6)
  return(output)
}


get_data_frame <- function(results){
  restriction <- c("NCMV", rep(" ", 2), "ACMV", rep(" ", 2))
  time <- c(paste0(seq(1, 3, by=1), " "), paste0(seq(1, 3, by=1), " "))
  
  MCAR <- c(results[[1]]$estimated_values1[7], results[[1]]$estimated_values2[7], results[[1]]$estimated_values3[7],
            results[[4]]$estimated_values1[7], results[[4]]$estimated_values2[7], results[[4]]$estimated_values3[7])
  MAR <- c(results[[2]]$estimated_values1[7], results[[2]]$estimated_values2[7], results[[2]]$estimated_values3[7],
           results[[5]]$estimated_values1[7], results[[5]]$estimated_values2[7], results[[5]]$estimated_values3[7])
  NMAR <- c(results[[3]]$estimated_values1[7], results[[3]]$estimated_values2[7], results[[3]]$estimated_values3[7],
           results[[6]]$estimated_values1[7], results[[6]]$estimated_values2[7], results[[6]]$estimated_values3[7])
  
  MCAR_sd <- c(sqrt(results[[1]]$var1[7]), sqrt(results[[1]]$var2[7]), sqrt(results[[1]]$var3[7]),
               sqrt(results[[4]]$var1[7]), sqrt(results[[4]]$var2[7]), sqrt(results[[4]]$var3[7]))
  MAR_sd <- c(sqrt(results[[2]]$var1[7]), sqrt(results[[2]]$var2[7]), sqrt(results[[2]]$var3[7]),
               sqrt(results[[5]]$var1[7]), sqrt(results[[5]]$var2[7]), sqrt(results[[5]]$var3[7]))
  NMAR_sd <- c(sqrt(results[[3]]$var1[7]), sqrt(results[[3]]$var2[7]), sqrt(results[[3]]$var3[7]),
              sqrt(results[[6]]$var1[7]), sqrt(results[[6]]$var2[7]), sqrt(results[[6]]$var3[7]))
  
  effect_df <- data.frame(restriction = restriction,
                          time = time,
                          MCAR = format(MCAR, digits = 2),
                          MAR = format(MAR, digits = 2),
                          NMAR = format(NMAR, digits = 1))
  
  sd_df <- data.frame(restriction = restriction,
                      time = time,
                      MCAR_sd = format(MCAR_sd, digits = 2),
                      MAR_sd = format(MAR_sd, digits = 2),
                      NMAR_sd = format(NMAR_sd, digits = 2))
  
  output <- data.frame(restriction = restriction,
                       time = time,
                       MCAR = stringr::str_c(effect_df[[3]], " (",sd_df[[3]], ")"),
                       MAR = stringr::str_c(effect_df[[4]], " (",sd_df[[4]], ")"),
                       NMAR = stringr::str_c(effect_df[[5]], " (",sd_df[[5]], ")"))
  return(output)
}


get_table <- function(data_frame){
  output <- data_frame |> kableExtra::kable(booktabs = TRUE, format = "latex")
  
  return(output)
}


save <- function(table, missing_rate){
  file_name <- paste0(missing_rate,"_table.txt")
  path <- here::here("RCT","04_report","output","table", file_name)
  writeLines(table, path)
}


main()



main_analysis_graph <- function(){
  #library(ggplot2)
  
  #read missing 
  missing_type <- "MCAR"
  restriction <- "NCMV"
  delta <- 0
  
  estimated_values <- read_results(missing_type, restriction, delta) 
  treatment_effects <- extract_effect(estimated_values)
  
  plot <- get_graph(treatment_effects, missing_type)
  

  return(plot)
}


read_results <- function(missing_type, restriction, delta){
  path1 <- here::here("RCT", "03_analyze", "output", restriction,
                      paste0(missing_type,"0.1","delta",delta,".obj" ))
  path2 <- here::here("RCT", "03_analyze", "output", restriction,
                      paste0(missing_type,"0.2","delta",delta,".obj" ))
  path3 <- here::here("RCT", "03_analyze", "output", restriction,
                      paste0(missing_type,"0.3","delta",delta,".obj" ))
  
  results1 <- readRDS(path1)
  results2 <- readRDS(path2)
  results3 <- readRDS(path1)
  
  output <- list(results1, results2, results3)
  return(output)
}


extract_effect <- function(results){
  effect1 <- c(results[1][[1]]$estimated_values1[7], results[1][[1]]$estimated_values2[7], results[1][[1]]$estimated_values3[7])
  effect2 <- c(results[2][[1]]$estimated_values1[7], results[2][[1]]$estimated_values2[7], results[2][[1]]$estimated_values3[7])
  effect3 <- c(results[3][[1]]$estimated_values1[7], results[3][[1]]$estimated_values2[7], results[3][[1]]$estimated_values3[7])
  output <- data.frame("estimated_values" = c(effect1, effect2, effect3),
                       "Time" = rep(seq(1,3, by= 1),3),
                       "missing_rate" = c(rep("0.1",3), rep("0.2",3), rep("0.3",3)))
  return(output)
}


get_graph <- function(treatment_effects, missing_type){
  real_parameter <- data.frame("estimated_values" = c(9, 6, 2),
                               "Time" = seq(1,3, by= 1),
                               "missing_rate" = rep("real", 3))
  treatment_effects <- treatment_effects |> dplyr::bind_rows(real_parameter)
  plot <- ggplot(data = treatment_effects, aes(x = Time, y = estimated_values, color = missing_rate)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("darkslategray2","deepskyblue2","dodgerblue4","coral2"))+
    scale_x_continuous(breaks = seq(0, 4, by = 1)) +
    scale_y_continuous(breaks = seq(0, 10, by = 2)) +
    labs(title = missing_type) +
    ylab("Treatment effect")+
    xlab("Time")
  
  return(plot)
}


plot <- main_analysis_graph()
plot
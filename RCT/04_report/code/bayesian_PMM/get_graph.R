get_graph <- function(){
  library(ggplot2)
  
  #read missing 
  missing_type <- "MCAR"
  
  #get graph
  estimated_values <- read_results(missing_type) 
  treatment_effects <- extract_effect(estimated_values)
  
  plot <- get_plot(treatment_effects, missing_type)
  
  
  #save plot
  # file_name <- paste0(missing_type, "delta", delta,".pdf")
  # path <- here::here("RCT", "04_report","output", restriction, file_name)
  # ggsave(filename = path, plot = plot, device = "pdf",width = 3, height = 5)
  
  return(plot)
}

read_results <- function(missing_type){
  path1 <- here::here("RCT", "03_analyze", "output", "bayesian_PMM",
                      paste0(missing_type,"0.1",".obj" ))
  path2 <- here::here("RCT", "03_analyze", "output", "bayesian_PMM",
                      paste0(missing_type,"0.2",".obj" ))
  path3 <- here::here("RCT", "03_analyze", "output", "bayesian_PMM",
                      paste0(missing_type,"0.3",".obj" ))
  
  results1 <- readRDS(path1)
  results2 <- readRDS(path2)
  results3 <- readRDS(path3)
  
  output <- list(results1, results2, results3)
  return(output)
}

extract_effect <- function(estimated_values){
  effect1 <- c(estimated_values[[1]]$beta1[2],estimated_values[[2]]$beta1[2],estimated_values[[3]]$beta1[2])
  effect2 <- c(estimated_values[[1]]$beta2[2],estimated_values[[2]]$beta2[2],estimated_values[[3]]$beta2[2])
  effect3 <- c(estimated_values[[1]]$beta3[2],estimated_values[[2]]$beta3[2],estimated_values[[3]]$beta3[2])
  output <- data.frame("estimated_values" = c(effect1, effect2, effect3),
                       "missing_rate" = paste0(rep(seq(0.1, 0.3, by= 0.1),3)),
                       "Time" = c(rep(1, 3), rep(2, 3), rep(3, 3)))
  return(output)
}


get_plot <- function(treatment_effects, missing_type){
  real_parameter <- data.frame("estimated_values" = c(9, 6, 2),
                               "Time" = seq(1,3, by= 1),
                               "missing_rate" = rep("real", 3))
  treatment_effects <- treatment_effects |> dplyr::bind_rows(real_parameter)
  plot <- ggplot(data = treatment_effects, aes(x = Time, y = estimated_values, color = missing_rate)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = c("darkslategray2","deepskyblue2","dodgerblue4","coral2"))+
    theme(panel.grid.minor.y = element_blank(),
          plot.title = element_text(vjust = -2, hjust = 0, size = 17))+
    scale_x_continuous(breaks = seq(0, 4, by = 1)) +
    scale_y_continuous(breaks = seq(-2, 10, by = 1), limits=c(-2,10)) +
    labs(title = missing_type) +
    ylab("Treatment effect")+
    xlab("Time")
}

plot <- get_graph()
plot

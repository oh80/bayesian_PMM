main <- function(){
  library(ggplot2)
  # read sample
  missing_rate <- 0.5
  missing_type <- "MCAR"
  file_name <- paste0(missing_type , "_", missing_rate, ".obj")
  path <- here::here("Liner_regression", "03_analyze","output", "bayesian_PMM", file_name)
  sample <- readRDS(path)
  
  # draw a graph
  beta <- 10
  missing_value <- 1
  plot <- draw_trace_plot(sample, beta, missing_value, missing_type)
  
  #save plot
  file_name <- paste0(missing_type, "trace_plot.pdf")
  path <- here::here("Liner_regression", "04_report", "output","bayesian_PMM",file_name)
  ggsave(filename = path, plot = plot, device = "pdf",  width = 8, height = 4)
}



draw_trace_plot <- function(sample, beta, missing_value, missing_type){
  sample1 <- sample[[missing_value]]
  sample1_beta <- sample1$beta_sample|> t()
  data <- data.frame("sample" = sample1_beta[,beta], "times" = 1:1000)
  output <- ggplot(data = data, mapping = aes(x = times, y = sample))+
    geom_line(colour = "seagreen")+
    labs(title = missing_type)
    
    
  return(output)
}


main()

main <- function(){
  library(ggplot2)
  # read sample
  missing_rate <- 0.5
  missing_type <- "MAR"
  file_name <- paste0(missing_type , "_", missing_rate, ".obj")
  path <- here::here("Liner_regression", "03_analyze","output", "bayesian_PMM", file_name)
  sample <- readRDS(path)
  
  # draw a graph
  plot <- draw_trace_plot(sample)
  return(plot)
}



draw_trace_plot <- function(sample){
  sample1 <- sample[[1]]
  sample1_beta <- sample1$beta_sample|> t()
  data <- data.frame("sample" = sample1_beta[,1], "times" = 1:1000)
  output <- ggplot(data = data, mapping = aes(x = times, y = sample))+
    geom_line()
  return(output)
}

plot <- main()

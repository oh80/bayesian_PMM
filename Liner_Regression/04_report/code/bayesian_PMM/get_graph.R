main <- function(){
  library("ggplot2")
  source("Liner_Regression/04_report/code/Bayesian_PMM/get_estimated_values.R")
  
  #get estimated values
  missing_type <- "NMAR"
  data <- store_estimated_values(missing_type)
  
  #get graph
  plot <- get_graph(missing_type, data)
  
  #save plot
  file_name <- paste0(missing_type, "estimated_values.pdf")
  path <- here::here("Liner_regression", "04_report", "output","bayesian_PMM",file_name)
  ggsave(filename = path, plot = plot, device = "pdf",  width = 8, height = 4)
}


store_estimated_values <- function(missing_type){
  estimated_values1 <- get_estimated_values(missing_rate = 0.1, missing_type = missing_type)
  estimated_values2 <- get_estimated_values(missing_rate = 0.25, missing_type = missing_type)
  estimated_values3 <- get_estimated_values(missing_rate = 0.5, missing_type = missing_type)
  
  output <- list()
  output[[1]] <- estimated_values1
  output[[2]] <- estimated_values2
  output[[3]] <- estimated_values3
  return(output)
}

  
get_graph <- function(missing_type ,data){
  x_axis <- c(rep(seq(1, 10, by = 1),4))
  estimated_values <- c(data[[1]] |> t(), data[[2]] |> t(), data[[3]] |> t(), 
                        seq(1, 10, by = 1))
  missing_rate <- c(rep("missing rate = 0.1", 10), rep("missing rate = 0.25", 10),
                    rep("missing rate = 0.5", 10), rep("real", 10))
  
  data <- data.frame("estimated_values" = estimated_values, 
                     "x_axis" = x_axis, 
                     "missing_rate" = missing_rate )
  
  plot <- ggplot(data = data, aes(x = x_axis, y = estimated_values, color = missing_rate)) +
    geom_point() +
    scale_x_continuous(breaks = seq(1, 10, by = 1)) +
    scale_y_continuous(breaks = seq(1, 10, by = 1)) +
    scale_color_manual(values = c("darkslategray2","deepskyblue2","dodgerblue4","coral2"  ))+
    geom_abline(intercept = 0,slope = 1,linetype=1,color= "coral2" ) +
    labs(title = missing_type) +
    ylab("beta")+
    xlab(" ")
  
  return(plot)
}


main()


main <- function(){
  library("ggplot2")
  
  #read data 
  missing_type <- "NMAR"
  data <- read_data(missing_type)
  
  #get graph
  plot <- get_graph(missing_type, data)
  
  #save plot
  file_name <- paste0(missing_type, "estimated_values.pdf")
  path <- here::here("Liner_regression", "04_report","output", "CCMV",file_name)
  ggsave(filename = path, plot = plot, device = "pdf",  width = 8, height = 4)
}


read_data <- function(missing_type){
  file_name1 <- paste0(missing_type,"_0.1.obj")
  file_name2 <- paste0(missing_type,"_0.25.obj")
  file_name3 <- paste0(missing_type,"_0.5.obj")
  
  data1 <- readRDS(here::here("Liner_regression", "03_analyze","output", "CCMV",file_name1))
  data2 <- readRDS(here::here("Liner_regression", "03_analyze","output", "CCMV",file_name2))
  data3 <- readRDS(here::here("Liner_regression", "03_analyze","output", "CCMV",file_name3))
  
  output <- list()
  output[[1]] <- data1
  output[[2]] <- data2
  output[[3]] <- data3
  return(output)
}


get_graph <- function(missing_type ,data){
  x_axis <- c(rep(seq(1, 10, by = 1),4))
  estimated_values <- c(data[[1]]$estimated_values[2:11], data[[2]]$estimated_values[2:11],
                        data[[3]]$estimated_values[2:11], 
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

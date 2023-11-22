main <- function(){
  library(ggplot2)
  
  #read data
  path <- here::here("RCT","01_data","data","data.obj")
  data <- readRDS(path)
  
  #get graph
  plot <- data |> extract_outcome() |> get_plot()
  
  #save
  file_name <- paste0("raw_data_outcome",".pdf")
  path <- here::here("RCT", "04_report","output","complete_data",file_name)
  ggsave(filename = path, plot = plot, device = "pdf",width = 5, height = 4)
  
}


extract_outcome <- function(data){
  mean <- data |> dplyr::group_by(Treatment) |>
    dplyr::summarise(Y1 = mean(Y1),
                     Y2 = mean(Y2),
                     Y3 = mean(Y3),
                     Y0 = mean(Y0)) 
 
  output <- mean |> tidyr::gather(Treatment) |> 
    dplyr::mutate("Time" = c(1,1,2,2,3,3,0,0)) |> 
    dplyr::mutate("group" = rep(c("Control","Treatment"),4))
  
  return(output)
}


get_plot <- function(data){
  plot <- ggplot(data = data,
                 mapping = aes(x = Time, y = value, colour = group )) +
    geom_point(size = 2.2) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = c("forestgreen","orange")) +
    theme(panel.grid.minor.y = element_blank(),
          panel.grid.minor.x  = element_blank())+
    scale_x_continuous(breaks = seq(0, 3, by = 1), limits=c(-0.3,3.3)) +
    scale_y_continuous(breaks = seq(-2, 11, by = 1), limits=c(-2,11)) +
    ylab("Outcome")+
    xlab("Time")
    
  return(plot)
}


main()


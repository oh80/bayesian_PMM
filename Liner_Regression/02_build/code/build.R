main <- function(){
  #read data
  path <- here::here("Liner_regression","01_data","data","row_data.obj")
  raw_data <- readRDS(path)
  
  #set params
  set.seed(428)
  missing_rate <- 0.5
  
  #MCAR
  missing_indicater <- get_index(raw_data, missing_rate)
  data_MCAR <- get_MCAR(raw_data, missing_indicater)
  return(data_MCAR)
  
  
}


get_index <- function(raw_data, missing_rate){
  N <- length(raw_data$Y)
  missing_indicater <- seq(0, 10, by = 1)
  prob <- c(1-missing_rate, rep(missing_rate/10, 10))
  output <- sample(missing_indicater, N, prob, replace = TRUE)
  return(output)
}


get_MCAR <- function(raw_data, missing_indicater){
  data <- raw_data |> dplyr::mutate("R" = missing_indicater)
  for (i in 1:1000) {
    if(data[i,]$R >= 1){
      data[i,data[i,]$R] = NaN
    }
  }
  output <- data
  return(output)
}

main()



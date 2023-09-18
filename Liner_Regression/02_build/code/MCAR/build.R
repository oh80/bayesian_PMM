main <- function(){
  #read data
  path <- here::here("Liner_regression","01_data","data","row_data.obj")
  raw_data <- readRDS(path)
  
  #set params
  set.seed(428)
  missing_rate <- 0.25
  
  #get missing data
  missing_indicater <- get_index(raw_data, missing_rate)
  missing_data <- get_missing_data(raw_data, missing_indicater)

  #save data
  save(missing_data, missing_rate)
}


get_index <- function(data, missing_rate){
  N <- length(data$Y)
  missing_indicater <- seq(0, 10, by = 1)
  prob <- c(1-missing_rate, rep(missing_rate/10, 10))
  output <- sample(missing_indicater, N, prob, replace = TRUE)
  return(output)
}


get_missing_data <- function(data, missing_indicater){
  N <- length(data$Y)
  data <- data |> dplyr::mutate("R" = missing_indicater)
  for (i in 1:N) {
    if(data[i,]$R >= 1){
      data[i,data[i,]$R] = NaN
    }
  }
  output <- data
  return(output)
}


save <- function(missing_data, missing_rate){
  file_name <- paste0("MCAR_",missing_rate,".obj")
  path <- here::here("Liner_regression", "02_build","data", file_name)
  saveRDS(missing_data, path)
}


main()


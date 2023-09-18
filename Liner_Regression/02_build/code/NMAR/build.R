main <- function(){
  #read data
  path <- here::here("Liner_regression","01_data","data","row_data.obj")
  raw_data <- readRDS(path)
  
  #set params
  set.seed(428)
  missing_rate <- 0.1
  
  data_and_indicater <- raw_data |> add_noise() |> get_index(missing_rate)
  missing_data <- get_missing_data(data_and_indicater)
  
  save(missing_data, missing_rate)
}


add_noise <- function(data){
  N <- length(data$Y)
  for (j in 1:10) {
    set.seed
    col_name <- paste0("x",j,"_noise")
    data <- data |> dplyr::mutate("{col_name}" := data[,j] + rnorm(n = N, mean = 0, sd = 1))
  }
output <- data
return(output)
}


get_index <- function(data, missing_rate){
  missing_num <- missing_rate * length(data$Y) /10
  output <- data |> dplyr::mutate("R" = 0)
  for (j in 1:10) {
    col_name <- paste0("x",j,"_noise")
    output <- output |> dplyr::arrange(R) |> dplyr::arrange(-.data[[col_name]])
    low_index <- 1 + missing_num*(j-1)
    high_index <- missing_num + missing_num*(j-1)
    output[low_index: high_index, ]$R <- j
  }
  return(output)
}


get_missing_data <- function(data_and_indicater){
  col_names <- c(paste0("x", seq(1, 10 ,by=1)), "Y", "R")
  output <- data_and_indicater |> dplyr::select(all_of(col_names))
  for (i in 1:1000) {
    if(output[i, ]$R >= 1){
      output[i,output[i,]$R] = NaN
    }
  }
  return(output)
}


save <- function(missing_data, missing_rate){
  file_name <- paste0("NMAR_",missing_rate,".obj")
  path <- here::here("Liner_regression", "02_build","data", file_name)
  saveRDS(missing_data, path)
}


main()

main <- function(){
  #read data
  path <- here::here("Liner_regression","01_data","data","row_data.obj")
  raw_data <- readRDS(path)
  
  #set params
  set.seed(428)
  missing_rate <- 0.5
  
  data_and_indicater <- raw_data |> add_col() |> get_index(missing_rate)
  missing_data <- get_missing_data(data_and_indicater)
  
  save(data_and_indicater, missing_rate)
}


add_col <- function(data){
  for (j in 1:10) {
    set.seed(j)
    epsilon <- rnorm(n = 1, mean = 0, sd = 1)
    col_name <- paste0("x",j,"_missing_predict")
    if(j <= 9){
      predict_col = j+1
      data <- data |> dplyr::mutate("{col_name}" := data[,predict_col] + epsilon )
    }else{
      data <- data |> dplyr::mutate("{col_name}" := data[,1] + epsilon )
    }
  }
  output <- data
  return(output)
}


get_index <- function(data, missing_rate){
  missing_num <- missing_rate * length(data$Y) /10
  output <- data |> dplyr::mutate("R" = 0)
  for (j in 1:10) {
    col_name <- paste0("x",j,"_missing_predict")
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
  file_name <- paste0("MAR_",missing_rate,".obj")
  path <- here::here("Liner_regression", "02_build","data", file_name)
  saveRDS(missing_data, path)
}


main()

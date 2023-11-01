main <- function(){
  path <- here::here("RCT","01_data","data","data.obj")
  data <- readRDS(path)
  
  set.seed(428)
  missing_rate <- 0.3
  pattern_prob <- get_prob(missing_rate)
  
  missing_data <- data |> get_missing_indicater(pattern_prob) |> 
    get_missing_data()
  
  save(missing_data, missing_rate)
}


get_prob <- function(missing_rate){
  prob0 <- missing_rate
  prob1 <- (1 - prob0) * missing_rate
  prob2 <- (1 - prob1 - prob0) * missing_rate
  prob3 <- (1 - prob2 - prob1 - prob0)
  
  output <- c(prob0, prob1, prob2, prob3)
  return(output)
}


get_missing_indicater <- function(data, pattern_prob){
  N <- length(data$X1)
  missing_indicater <- seq(0, 3, by = 1)
  R <- sample(missing_indicater, N, pattern_prob, replace = TRUE)
  
  output <- data |> dplyr::mutate("R" = R)
  return(output)
}


get_missing_data <- function(data){
  N <- length(data$X1)
  for (i in 1:N) {
    if(data$R[i] == 0){
      data$Y1[i] = data$Y2[i] = data$Y3[i] = NaN
    }
    if(data$R[i] == 1){
      data$Y2[i] = data$Y3[i] = NaN
    }
    if(data$R[i] == 2){
      data$Y3[i] = NaN
    }
  }
  output <- data
  return(output)
}


save <- function(missing_data, missing_rate){
  file_name <- paste0("MCAR_", missing_rate, ".obj")
  path <- here::here("RCT","02_build","data","MCAR",file_name)
  saveRDS(missing_data, path)
}


main()

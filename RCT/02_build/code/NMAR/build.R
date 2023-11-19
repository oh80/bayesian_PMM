main <- function(){
  path <- here::here("RCT","01_data","data","data.obj")
  data <- readRDS(path)
  
  set.seed(428)
  missing_rate <- 0.2
  missing_data <- get_missing_indicater(data, missing_rate) |> 
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


get_missing_indicater <- function(data, missing_rate){
  N <- length(data$X1)
  effect_for_missing <- c(2,2,2,2)
  epsilon <- rnorm(N , mean = 0, sd = 2)
  
  missing_predict <- rep(0, N)
  for (i in 1:N) {
    values <- c(data$Y0[i], data$Y1[i], data$Y2[i], data$Y3[i])
    missing_predict[i] <- sum(values * effect_for_missing) + epsilon[i]
  }
  
  data <- data |> dplyr::mutate("missing_predict" = missing_predict) |> 
    dplyr::arrange(-missing_predict)
  
  pattern_prob <- get_prob(missing_rate)
  missing_num <- N*pattern_prob
  
  data <- data |> dplyr::mutate("R" = rep(3,N))
  pattern_devid <- c(1,missing_num[1], missing_num[1]+1, missing_num[1]+missing_num[2],
                     missing_num[1]+missing_num[2]+1, missing_num[1]+missing_num[2]+missing_num[3])
  data$R[pattern_devid[1]: pattern_devid[2]] <- 0
  data$R[pattern_devid[3]: pattern_devid[4]] <- 1
  data$R[pattern_devid[5]: pattern_devid[6]] <- 2
  
  output <- data |> dplyr::select(-missing_predict)
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
  file_name <- paste0("NMAR_", missing_rate, ".obj")
  path <- here::here("RCT","02_build","data","NMAR",file_name)
  saveRDS(missing_data, path)
}


main()
main <- function(){
  set.seed(428)
  sample_size <- 1000
  
  baseline_data <- generate_baseline(sample_size)
  data <- baseline_data |> assign_group() |> 
    generate_follow_up()
  
  save(data)
}

generate_baseline <- function(sample_size){
  control_varable_mean <- rep(0, 5)
  control_varable_cov <- matrix(c(2.0, 0.5, 0.5, 0.5, 0.5,
                                  0.5, 2.0, 0.5, 0.5, 0.5,
                                  0.5, 0.5, 2.0, 0.5, 0.5,
                                  0.5, 0.5, 0.5, 2.0, 0.5,
                                  0.5, 0.5, 0.5, 0.5, 2.0), ncol = 5)
  control_varable_effect <- c(1, 2, 3, 4, 5)
  
  control_varable <- mvtnorm::rmvnorm(n = sample_size, control_varable_mean, control_varable_cov )
  
  baseline_data <- rep(0, sample_size)
  for (i in 1:sample_size) {
    baseline_data[i] <- sum(control_varable[i,] * control_varable_effect)
  }
  
  output <- control_varable |> as.data.frame() |> 
    magrittr::set_colnames(paste0("X", seq(1, 5, by = 1))) |> 
    dplyr::mutate("Y0" = baseline_data)
  
  return(output)
}


assign_group <- function(baseline_data){
  N <- length(baseline_data$X1)
  group <- c(rep(0, N/2), rep(1, N/2))
  output <- baseline_data |> dplyr::mutate("Treatment" = group)
  
  return(output)
}


generate_follow_up <- function(baseline_data){
  N <- length(baseline_data$X1)
  #generate first follow up
  data <- baseline_data |> dplyr::mutate("Y1" = rep(0, N))
  for(i in 1:N){
    if(data$Treatment[i] == 1){
      epsilon <- rnorm(1, 10, 2)
      data$Y1[i] = data$Y0[i] + epsilon
    }
    if(data$Treatment[i] == 0){
      epsilon <- rnorm(1, 1, 2)
      data$Y1[i] = data$Y0[i] + epsilon
    }
  }
  #generate second follow up
  data <- data |> dplyr::mutate("Y2" = rep(0, N))
  for(i in 1:N){
    if(data$Treatment[i] == 1){
      epsilon <- rnorm(1, -2, 2)
      data$Y2[i] = data$Y1[i] + epsilon
    }
    if(data$Treatment[i] == 0){
      epsilon <- rnorm(1, 1, 2)
      data$Y2[i] = data$Y1[i] + epsilon
    }
  }
  #generate third follow up
  data <- data |> dplyr::mutate("Y3" = rep(0, N))
  for(i in 1:N){
    if(data$Treatment[i] == 1){
      epsilon <- rnorm(1, -3, 2)
      data$Y3[i] = data$Y2[i] + epsilon
    }
    if(data$Treatment[i] == 0){
      epsilon <- rnorm(1, 1, 2)
      data$Y3[i] = data$Y2[i] + epsilon
    }
  }
  return(data)
}


save <- function(data){
  path <- here::here("RCT","01_data","data","data.obj")
  saveRDS(data, path)
}


 main()

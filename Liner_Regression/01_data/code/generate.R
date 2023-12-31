main <- function(){
  #set params
  set.seed(428)
  beta <- seq(1,10,by=1)
  sigma <- 0.5
  N = 1000
  
  #generate params at random
  mu_x <- get_mu()
  cov_x <- get_cov()
  
  #generate data
  X <- generate_x(N, mu_x, cov_x)
  Y <- generate_y(X, beta, sigma)
  
  #save_data
  save(X,Y)
}


get_mu <- function(){
  mu_0 <- rep(0, each = 10)
  output <- mvtnorm::rmvnorm(n = 1, mean = mu_0)
  
  return(output)
}


get_cov <- function(){
  W <- rWishart(n = 1, df =10, Sigma = diag(10)) |> matrix(ncol = 10)
  output <- rWishart(n = 1, df =10, Sigma = W)|> matrix(ncol = 10)
  
  return(output)
}


generate_x <- function(N, mu_x, cov_x){
  col_name <- paste0("x",seq(1,10,by=1))
  output <- mvtnorm::rmvnorm(n = N, mean = mu_x, sigma = cov_x) |>
    scale() |> as.data.frame() |> 
    magrittr::set_colnames(col_name)
  
  return(output)
}


generate_y <- function(X, beta, sigma){
  mu_y <- rep(0, length(X$x1))
  for (i in 1:length(mu_y)) {
    mu_y[i] <- sum(X[i,] * beta)
  }
  
  output <- rep(0,length(mu_y)) |>
    as.data.frame() |>
    magrittr::set_colnames("Y")
  for (i in 1:length(mu_y)) {
    output[i,1] <- rnorm(n = 1, mean = mu_y[i], sd = sigma)
  }
  
  return(output)
}


save <- function(X,Y){
  merge_data <- X |> dplyr::bind_cols(Y)
  path <- here::here("Liner_regression","01_data","data","row_data.obj")
  saveRDS(merge_data, path)
}


main()


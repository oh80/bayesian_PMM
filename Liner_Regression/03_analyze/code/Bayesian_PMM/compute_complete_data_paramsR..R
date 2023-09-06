get_OLS_extimater <- function(data){
  complete_data <- data |> dplyr::filter(R == 0) |> dplyr::select(-R)
  X <- complete_data |> dplyr::select(-Y)
  Y <- complete_data |> dplyr::select(Y)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  X_t <- t(X)
  X_square <- X_t %*% X
  
  output <- solve(X_square) %*% X_t %*% Y
  return(output)
}


get_Sigma0 <- function(data, g, sigma){
  complete_data <- data |> dplyr::filter(R == 0) |> dplyr::select(-R)
  X <- complete_data |> dplyr::select(-Y)
  
  X <- as.matrix(X)
  X_t <- t(X)
  X_square <- X_t %*% X
  
  output <- g * sigma * solve(X_square) 
  return(output)
}



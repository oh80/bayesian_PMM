data <- readRDS(here::here("Liner_regression", "03_analyze","output", "CCMV","MCAR_0.5.obj"))

data$estimation
data$var

data2 <- readRDS(here::here("Liner_regression", "03_analyze","output", "CCMV","MCAR_0.1.obj"))

data2$estimation
data2$var


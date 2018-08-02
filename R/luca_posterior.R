#' @importFrom mvtnorm dmvnorm 

#mu_true=[1.6 1.4]';
#Marglike_true=1;
lposterior_1 <- function(x){
  mu_1 <- c(-10, -10)
  sigma_1 <- matrix(c(2, rep(0.6, 2), 1), ncol = 2, nrow = 2)
  mu_2 <- c(0, 16)
  sigma_2 <- matrix(c(2, rep(-0.4, 2), 2), ncol = 2, nrow = 2)
  mu_3 <- c(13, 8)
  sigma_3 <- matrix(c(2, rep(0.8, 2), 2), ncol = 2, nrow = 2)
  mu_4 <- c(-9, 7)
  sigma_4 <- matrix(c(3, rep(0, 2), 5), ncol = 2, nrow = 2)
  mu_5 <- c(14, -14)
  sigma_5 <- matrix(c(2, rep(-0.1, 2), 2), ncol = 2, nrow = 2)
  f_1  <- dmvnorm(x, mean = mu_1 , sigma = sigma_1) 
  f_2  <- dmvnorm(x, mean = mu_2 , sigma = sigma_2) 
  f_3  <- dmvnorm(x, mean = mu_3 , sigma = sigma_3) 
  f_4  <- dmvnorm(x, mean = mu_4 , sigma = sigma_4) 
  f_5  <- dmvnorm(x, mean = mu_5 , sigma = sigma_5) 
  f <- f_1 + f_2 + f_3 + f_4 + f_5
  # print(paste("f_1 = " ,  f_1))
  # print(paste("f_2 = " ,  f_2))
  # print(paste("f_3 = " ,  f_3))
  # print(paste("f_4 = " ,  f_4))
  # print(paste("f_5 = " ,  f_5))
  log(f) - log(5)
}
  
# mu_true=[0 16]';
# Marglike_true=1;
# mu2=[0 16];
# SIGMA2 = [3 0;0 3];
lposterior_2 <- function(x){
  mu <- c(0, 16)
  sigma <- matrix(c(3, rep(0, 2), 3), ncol = 2, nrow = 2)
  dmvnorm(x, mean = mu , sigma = sigma, log = TRUE) 
}


# mu_true=[2.5 8];
# Marglike_true=1;
# mu1=[5 0];
# SIGMA1 = [2 0.6; 0.6 1];
# mu2=[0 16];
# SIGMA2 = [3 0;0 3];
# f1=1/2*mvnpdf(x,mu1,SIGMA1);
# f2=1/2*mvnpdf(x,mu2,SIGMA2);
# f=f1+f2;
lposterior_3 <- function(x){
  mu_1 <- c(5, 0)
  sigma_1 <- matrix(c(2, rep(0.6, 2), 1), ncol = 2, nrow = 2)
  mu_2 <- c(0, 16)
  sigma_2 <- matrix(c(3, rep(0, 2), 3), ncol = 2, nrow = 2)
  f_1  <- 1/2 * dmvnorm(x, mean = mu_1 , sigma = sigma_1) 
  f_2  <- 1/2 * dmvnorm(x, mean = mu_2 , sigma = sigma_2) 
  f <- f_1 + f_2 
  log(f)
}

# mu_true=[0 16 5 -5]';
# Marglike_true=1;
# mu2=[0 16 5 -5];
# sig=4;
# SIGMA2 = sig*eye(DIM);
# f=mvnpdf(x,mu2,SIGMA2);
lposterior_4 <- function(x){
  mu <- c(0, 16, 5, -5)
  sigma <- diag(4, 4)
  dmvnorm(x, mean = mu , sigma = sigma, log = TRUE) 
}

# mu <- true=5*ones(DIM,1);
# Marglike <- true=1;
# mu2=5*ones(1,DIM);
# sig=4;
# SIGMA2 = sig*eye(DIM);
# f=mvnpdf(x,mu2,SIGMA2);
lposterior_5 <- function(x){
  mu <- rep(5, 10)
  sigma <- diag(4, 10)
  dmvnorm(x, mean = mu , sigma = sigma, log = TRUE) 
}


lposterior_6 <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  logtarget <-  -1/32 * (4 - 10 * x1 - x2^2)^2 - x1^2/50 - x2^2/50
}


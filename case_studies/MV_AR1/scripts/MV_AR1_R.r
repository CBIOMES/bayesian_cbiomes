library(MASS)     #package for the multi-variate normal distribution
library(rstan)
options(mc.cores = parallel::detectCores())

p     <- 3

SIGMA  <- diag(runif(p))   #generate random numbers

eig    <- runif(p,0,1)
Z      <- matrix(ncol=p, rnorm(p^2))
decomp <- qr(Z)
Q      <- qr.Q(decomp)
R      <- qr.R(decomp)
d      <- diag(R)
ph     <- d / abs(d)
O      <- Q %*% diag(ph)
PHI    <- t(O) %*% diag(eig) %*% O

T  <- 200
y0 <- rnorm(p)
Y     <- matrix(NA,p,T)
Y[,1] <- y0

for(t in 2:T){
    Y[,t] <- PHI%*%Y[,t-1] + mvrnorm(1,rep(0,3),SIGMA)}

options(repr.plot.width=6, repr.plot.height=4)
matplot(t(Y),type='l')

mod_code <- "data {
	int T;         //length of time series
	int p;         //number of variables
	matrix[p,T] Y; //matrix of observations; variables are rows; time is columns
}
parameters{
	matrix[p,p] PHI;     //dynamics matrix
	vector<lower=1E-15>[p] sigma;     //variances of stochastic forcing
	vector[p] init;      //mean of initial conditions as parameter vector
}
model{
	Y[,1] ~ normal(init, sigma);           //distribution of the initial conditions
	for(i in 2:T){
		Y[,i] ~ normal(PHI*Y[,i-1],sigma); //conditional predictive distribution
	}
}"

mod <- stan_model(model_code=mod_code)

data <- list(p=p,T=T,Y=Y)

mcmc <- sampling(mod,data=data,iter=2000,warmup=1000,open_progress=TRUE)

mcmc

A      <- matrix(runif(p^2)*2-1, ncol=p)   #generate random numbers
SIGMA2 <- t(A) %*% A   
SIGMA2

T      <- 200
y20    <- rnorm(p)
Y2     <- matrix(NA,p,T)
Y2[,1] <- y20

for(t in 2:T){
    Y2[,t] <- PHI%*%Y2[,t-1] + mvrnorm(1,rep(0,3),SIGMA2)}

options(repr.plot.width=6, repr.plot.height=4)
matplot(t(Y2),type='l')

data2 <- list(p=p,T=T,Y=Y2)

mod_code_cov <- "data {
	int T;         //length of time series
	int p;         //number of variables
	matrix[p,T] Y; //matrix of observations; variables are rows; time is columns
}
parameters{
	matrix[p,p] PHI;     //dynamics matrix
	cov_matrix[p] SIGMA; //co-variance matrix of stochastic forcing
	vector[p] init;      //mean of initial conditions as parameter vector
}
model{
	Y[,1] ~ multi_normal(init, SIGMA);           //distribution of the initial conditions
	for(i in 2:T){
		Y[,i] ~ multi_normal(PHI*Y[,i-1],SIGMA); //conditional predictive distribution
	}
}"

mod_cov <- stan_model(model_code=mod_code_cov)

mcmc_cov <- sampling(mod_cov,data=data2,chains=4,open_progress=TRUE)

mcmc_cov

mod_code_D_struc <- "data {
	int T;         //length of time series
	int p;         //number of variables
	matrix[p,T] Y; //matrix of observations; variables are rows; time is columns
}
parameters{
	matrix[p,p] PHI;               //dynamics matrix
	vector<lower=1E-15>[p] sigma;  //variances of stochastic forcing
	vector[p] init;                //mean of initial conditions as parameter vector
}
model{
    PHI[1,3] ~ normal(0,1E-3);
    PHI[3,1] ~ normal(0,1E-3);

	Y[,1] ~ normal(init, sigma);           //distribution of the initial conditions
	for(i in 2:T){
		Y[,i] ~ normal(PHI*Y[,i-1],sigma); //conditional predictive distribution
	}
}"

mod_struc <- stan_model(model_code=mod_code_D_struc)

mcmc_struc <- sampling(mod_struc,data=data,chains=4,open_progress=TRUE)

mcmc_struc

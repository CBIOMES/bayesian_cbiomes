library(deSolve)
library(rstan)
options(mc.cores = parallel::detectCores())

##########################################################

theta <- list(gamma =0.25, lambda=0.1)
x <- c(P=1)
T  <- 30
dt <- 1
nt <- T/dt
t  <- seq(0,T,length.out=nt)

###########################################################

dxdt <- function(t,x,theta){
    with(as.list(c(x,theta)),{
        dP <- gamma*P - lambda*P*P
        list(c(dP)) })}
		
###########################################################

x <- as.data.frame(ode(y=x, times=t, func=dxdt, parms=theta))

###########################################################

plot(t,x$P,type='l')

###########################################################

t_obs_ind <- c(1,3,4,6,9,13,16,24)
obs       <- x$P[t_obs_ind] + rnorm(length(t_obs_ind),sd=0.1)
t_obs     <- x$t[t_obs_ind] 
plot(t_obs,obs)

###########################################################

stancode <- "functions {
   real[] P_growth(real t,       // time
                   real[] x,      // state
                   real[] theta, // parameters
                   real[] x_r,   // environmental data
                   int[] x_i){
    real gamma  = theta[1];
    real lambda = theta[2];

    real growth = gamma*x[1];
    real loss   = lambda*x[1]*x[1];
    
    return {growth - loss};
  }
}
data {
    int<lower = 0> N;           // num obs
    real t_obs[N];              // obs times
    real<lower = 0> y[N];       // observed variable at measurement times
    real sigma;
}
parameters {
    real<lower=0,upper=1> theta[2];      // parameters
    real<lower=0> x0[1];
}
transformed parameters {
    real x[N,1] = integrate_ode_rk45(P_growth, x0, -1, t_obs, theta,
                                      rep_array(0.0, 0), rep_array(0, 0),
                                      1e-6, 1e-5, 1e3);
}
model {
    //theta[1] ~ normal(0.1, 2);
    //theta[2] ~ normal(0.1, 2);
    x0       ~ normal(1.0, 10);
    y[1:N]   ~ normal(x[1:N,1], sigma); // obs
}"

########################################################

mod <- stan_model(model_code=stancode)

########################################################

data <- list(N=length(t_obs),
             t_obs=t_obs,
             y=obs,
             sigma=0.1)
			 
########################################################

mcmc <- sampling(mod,data=data,iter=2000,chains=4,open_progress=TRUE)

########################################################

mcmc

########################################################

post <- extract(mcmc)

########################################################

names(post)

########################################################

dim(post$theta)
dim(post$x0)

########################################################

par(mfrow=c(1,3))
hist(post$theta[,1],xlab='gamma',main='')
hist(post$theta[,2],xlab='lambda',main='')
hist(post$x0,xlab='x0',main='')

########################################################

plot(post$theta[,1],post$theta[,2])

########################################################

theta_sin <- list(gamma =0.25,
                  lambda=0.1,
                  period=365,
                  omega =2*pi/365,
                      dt=1)
x_sin <- c(P   =2.5,
           time=0)
T_sin  <- 365*4
dt_sin <- 1
nt_sin <- T_sin/dt_sin

t_sin  <- seq(0,T_sin,length.out=nt_sin)

#######################################################

dxdt_sin <- function(t,x,theta){
    with(as.list(c(x,theta)),{
        dP      <- gamma*(1+sin(omega*time))*P - lambda*P*P
        delta_t <- dt
        list(c(dP, delta_t))
    })
}

#######################################################

x_sin <- as.data.frame(ode(y=x_sin, times=t_sin, func=dxdt_sin, parms=theta_sin))

#######################################################

plot(x_sin$time,x_sin$P,type='l')

#######################################################

t_obs_ind_sin      <- sort(sample(1:length(t_sin),50))
obs_sin            <- x_sin$P[t_obs_ind_sin] + rnorm(length(t_obs_ind_sin),sd=0.5)
obs_sin[obs_sin<0] <- 0
t_obs_sin          <- x_sin$time[t_obs_ind_sin]

#######################################################

plot(t_obs_sin,obs_sin,type='l')

#######################################################

stancode_sin <- "functions {
   real[] P_growth(real t,       // time
                   real[] x,      // state
                   real[] theta, // parameters
                   real[] x_r,   // environmental data
                   int[] x_i){
    real gamma  = theta[1];
    real lambda = theta[2];
	//real a      = theta[3];

    real growth = gamma*x[1] + gamma*sin(2*pi()*(1.0/365.0)*t)*x[1];
    real loss   = lambda*x[1]*x[1];
    
    return {growth - loss};
  }
}
data {
    int<lower = 0> N;           // num obs
    real<lower = 0> t_obs[N];              // obs times
    real<lower = 0> y[N];       // observed variable at measurement times
    //real<lower = 0> sigma;
}
parameters {
    real<lower=0> theta[2];      // parameters
    real<lower=0> x0[1];
	real<lower=1E-15> sigma;
}
transformed parameters {
    real<lower=0> x[N,1] = integrate_ode_rk45(P_growth, x0, 1, t_obs, theta,
                                      rep_array(0.0, 0), rep_array(0, 0),
                                      1e-6, 1e-6, 1e5);
}
model {
    theta[1] ~ normal(0.1, 1);
    theta[2] ~ normal(0.1, 1);
	//theta[3] ~ normal(0.01,1);
    x0       ~ normal(1.0, 10);
    y[1:N]   ~ normal(x[1:N,1], sigma); // obs
}"

#####################################################

mod_sin <- stan_model(model_code=stancode_sin)

#####################################################

data_sin <- list(N=length(t_obs_sin),
             t_obs=round(t_obs_sin),
             y=obs_sin)
			 
#####################################################

mcmc_sin <- sampling(mod_sin,data=data_sin,iter=2000,chains=4,open_progress=TRUE)

#####################################################

mcmc_sin

#####################################################

post_sin <- extract(mcmc_sin)

#####################################################

par(mfrow=c(1,3))
hist(post_sin$theta[,1],xlab='gamma',main='')
hist(post_sin$theta[,2],xlab='lambda',main='')
hist(post_sin$x0,xlab='x0',main='')

#####################################################

pairs(post_sin$theta)

#####################################################

mu  <- colMeans(post_sin$x[,,1])
std <- apply(post_sin$x[,,1],2,sd)

#####################################################

plot(t_obs_sin,colMeans(post_sin$x[,,1]),type='l')
points(t_obs_sin,obs_sin)


library(rstan)
options(mc.cores = parallel::detectCores())

N     <- 100
beta1 <- 1.5
beta0 <- 2
sigma <- 1.25

x <- rnorm(N)
y <- beta0 + beta1*x + rnorm(N,mean=0,sd=sigma)

options(repr.plot.width=5, repr.plot.height=5)
plot(x,y)

stancode <- "
data {
    int       N;
    vector[N] x;
    vector[N] y;
}
parameters {
    real              beta0;
    real              beta1;
    real<lower=1E-15> sigma;
}
model{
    // Priors
    beta0 ~ normal(0,100);
    beta1 ~ normal(0,100);
    
    // Likelihood
    y ~ normal(beta0 + beta1*x, sigma);
}"

mod <- stan_model(model_code=stancode)

data <- list(N=N,x=x,y=y)

mcmc <- sampling(mod, data=data, chains=4, iter=2000, open_progress=TRUE)

mcmc

post <- extract(mcmc)

names(post)

options(repr.plot.width=6, repr.plot.height=5)
par(mfrow=c(3,1),mar=c(4,4,1,1))
plot(post$beta0,ylab='beta0')
plot(post$beta1,ylab='beta1')
plot(post$sigma,ylab='sigma')

options(repr.plot.width=6, repr.plot.height=2.5)
par(mfrow=c(1,3))
hist(post$beta0,main='beta0',xlab='')
    abline(v=beta0,lwd=2)
hist(post$beta1,main='beta1',xlab='')
    abline(v=beta1,lwd=2)
hist(post$sigma,main='sigma',xlab='')
    abline(v=sigma,lwd=2)

options(repr.plot.width=4, repr.plot.height=4)
plot(x,y)
abline(mean(post$beta0), mean(post$beta1),lwd=2)
for(i in sample(1:length(post$beta0),20)){
    abline(post$beta0[i],post$beta1[i],lty=2,col=adjustcolor('black',alpha.f=0.3))
}

stancode_prior <- "
data {
    int       N;
    vector[N] x;
    vector[N] y;
}
parameters {
    real              beta0;
    real              beta1;
    real<lower=1E-15> sigma;
}
model{
    // Priors
    beta0 ~ normal(0,100);
    beta1 ~ normal(2.5,0.1);
    
    // Likelihood
    y ~ normal(beta0 + beta1*x, sigma);
}"

xin <- seq(0,3,0.01)
plot(xin,dnorm(xin,mean=2.5,sd=0.1),type='l',xlab='beta1',ylab='Density')

mod_prior <- stan_model(model_code=stancode_prior)

mcmc_prior <- sampling(mod_prior, data=data, chains=4, iter=2000, open_progress=TRUE)

mcmc_prior

post_prior <- extract(mcmc_prior)

options(repr.plot.width=6, repr.plot.height=2.5)
par(mfrow=c(1,3))
hist(post_prior$beta0,main='beta0',xlab='',freq=FALSE)
    abline(v=beta0,lwd=2)
hist(post_prior$beta1,main='beta1',xlab='',xlim=c(0.5,3),freq=FALSE)
    abline(v=beta1,lwd=2)
    lines(xin,dnorm(xin,mean=2.5,sd=0.1))
hist(post_prior$sigma,main='sigma',xlab='',freq=FALSE)
    abline(v=sigma,lwd=2)

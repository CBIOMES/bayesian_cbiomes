import pystan
import numpy as np

###################################################################

N = 100
beta1 = 1.5
beta0 = 2.0
sigma = 1.25

x = np.random.normal(size=N)
y = beta0 + beta1 * x + np.random.normal(size=N, scale=sigma)

###################################################################

import matplotlib.pyplot as plt

fig, ax = plt.subplots(subplot_kw={'aspect':1})
ax.scatter(x,y)

###################################################################

stan_code = '''data {
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
}
'''

###################################################################

mod = pystan.StanModel(model_code=stan_code)

###################################################################

data = {'N':N, 'x':x, 'y':y}

###################################################################

mcmc = mod.sampling(data=data, chains=4, iter=2000)

###################################################################

print(mcmc)

###################################################################

print(mcmc.flatnames)

###################################################################

fig, axs = plt.subplots(nrows=3, sharex=True, figsize=(12,10))
for ax,name in zip(axs,mcmc.flatnames):
    ax.plot(mcmc[name], linestyle='none', marker='.')
    ax.set(ylabel=name, xlabel='index')
	
###################################################################

true_value = {'beta0':beta0, 'beta1':beta1, 'sigma':sigma}

fig, axs = plt.subplots(ncols=3, figsize=(12,6))
for ax,name in zip(axs,mcmc.flatnames):
    ax.hist(mcmc[name])
    ax.axvline(true_value[name], color='black')
    ax.set(ylabel='frequency', title=name)
	
###################################################################

def abline(ax, intercept, slope, x=None, **kwargs):
    if x is None:
        x = np.array(ax.get_xlim())
    ax.plot(x, intercept+slope*x, **kwargs)

fig, ax = plt.subplots()
ax.scatter(x,y)
xlim = np.array(ax.get_xlim())
abline(ax, intercept=np.mean(mcmc['beta0']), slope=np.mean(mcmc['beta1']), x=xlim, color='black', linewidth=2)
for i in range(20):
    abline(ax, intercept=mcmc['beta0'][i], slope=mcmc['beta1'][i], x=xlim, color='black', alpha=0.1)
	
##################################################################

stan_code_prior = '''data {
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
}
'''

##################################################################

import scipy.stats as stats

xin = np.linspace(0,3,200)

fig, ax = plt.subplots()
ax.plot(xin, stats.norm.pdf(xin, loc=2.5, scale=0.1))
ax.set(xlabel='beta1', ylabel='density')
None

##################################################################

mod_prior = pystan.StanModel(model_code=stan_code_prior)

##################################################################

mcmc_prior = mod_prior.sampling(data=data, chains=4, iter=2000)

##################################################################

print(mcmc_prior)

##################################################################

fig, axs = plt.subplots(ncols=3, figsize=(12,6))
for ax,name in zip(axs,('beta0','beta1','sigma')):
    ax.hist(mcmc_prior[name], density=True)
    ax.axvline(true_value[name], color='black')
    ax.set(ylabel='density', title=name)
axs[1].plot(xin, stats.norm.pdf(xin, loc=2.5, scale=0.1))





#!/usr/bin/env python
# coding: utf-8

# ## Introduction to Stan via linear regression
# 
# In this notebook we will use a simple example of single-variable linear regression to demonstrate the basic structure of a Stan program, how it is called, and how to plot the results Stan yields.
# 
# First we load Stan and tell Stan to look for and use multiple processors, if available.

# In[1]:


import pystan
import numpy as np


# The code below generates a synthetic dataset of length 100, with slope of 1.5, instercept of 2, and a noise error standard deviation of 1.25.
# $x$ is the independent variable and $y$ is the dependent variable.

# In[2]:


N = 100
beta1 = 1.5
beta0 = 2.0
sigma = 1.25

x = np.random.normal(size=N)
y = beta0 + beta1 * x + np.random.normal(size=N, scale=sigma)


# We plot the data

# In[19]:


import matplotlib.pyplot as plt

fig, ax = plt.subplots(subplot_kw={'aspect':1})
ax.scatter(x,y)


# The code below is our first Stan program. 
# There are several ways to pass Stan code to the Stan program.
# Here we are using character strings to easily integrate with a standalone notebook.
# 
# Stan programs are structured as templated blocks. 
# 
# #### `data` block
# We pass the data and all inputs in the `data` block. 
# In this case Stan needs the vectors $x$ and $y$ and the length of those vectors $N$.
# 
# #### `parameter` block
# We tell Stan which unknown parameters we seek to estimate via the `parameters` block.
# We will estimate the intercept `beta0`, the slope `beta1` and the error standard devation `sigma`.
# We give `sigma` a lower bound of a very small positive number to tell Stan not to try fitting with a standard deviation less than or equal to zero because it doesn't make mathemtical sense and Stan will get angry.
# 
# #### `model` block
# The `model` block contains two components (despite it only being one block): the prior and likelihood.
# This is where we tell Stan about the probability distributions being used in the problem.
# You can tell the prior and likelihood apart in that the prior will not use the data, while the likelihood will. 
# The likelihood will describe the distribution of the data *conditional on the parameters*.
# Here we assign normal distributions with mean zero and standard deviation 100 to `beta0` and `beta1`.
# We assign a normal distribution to the data with mean `beta0 + beta1*x` and standard deviation `sigma`.

# In[4]:


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


# We compile the Stan code below.

# In[5]:


mod = pystan.StanModel(model_code=stan_code)


# Organize the data for PyStan

# In[6]:


data = {'N':N, 'x':x, 'y':y}


# Perform mcmc

# In[7]:


mcmc = mod.sampling(data=data, chains=4, iter=2000)


# We look at the summary table for the mcmc results.

# In[8]:


print(mcmc)


# (See slides for a description of `n_eff` and `Rhat` diagnostics)

# In[9]:


print(mcmc.flatnames)


# The order of the mcmc samples represents the order in which they were sampled. 
# This means you can plot the mcmc samples as a sequence to see the trajectory of the sampler.
# 
# You can compare the appearance of this mcmc trajectory to that of the random walk MH algorithm. 
# In this case the samples are uncorrelated so each sample contributes much more information that the strongly correlated sample characteristic of RWMH.

# In[10]:


fig, axs = plt.subplots(nrows=3, sharex=True, figsize=(12,10))
for ax,name in zip(axs,mcmc.flatnames):
    ax.plot(mcmc[name], linestyle='none', marker='.')
    ax.set(ylabel=name, xlabel='index')


# Below we ignore the sequence and plot the distribution as a histogram.
# These plots represent the marginal distribution of each parameter after integrating (via mcmc) over the distribution of the other variables. 

# In[11]:


true_value = {'beta0':beta0, 'beta1':beta1, 'sigma':sigma}

fig, axs = plt.subplots(ncols=3, figsize=(12,6))
for ax,name in zip(axs,mcmc.flatnames):
    ax.hist(mcmc[name])
    ax.axvline(true_value[name], color='black')
    ax.set(ylabel='frequency', title=name)


# Below we randomly sample the mcmc samples to plot some lines.
# Since the order of the samples represents the mcmc steps, we want to pair the same sample numbers in order to get a proper sample from the joint posterior.

# In[12]:


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


# Previously we gave very wide normal distributions which is effectively *un*informative in the sense that the solution is indistinguishable had we given the uniform priors.
# 
# EXERCISE IDEA: try imposing a uniform prior in the above and compare the fits.
# 
# Below we demonstrate the effect of *informative* prior distributions. 
# We give the prior distribution on the slope a normal distribution with mean 2.5 and standard deviation of 0.1.
# 
# For example, this prior may have come from analysis of previous data. 
# The Bayesian solution has a nice property that we would get the same answer when we use a previous posterior as a prior, vs. if we pooled all the data and analyzed them together.

# In[13]:


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


# We plot the prior.

# In[14]:


import scipy.stats as stats

xin = np.linspace(0,3,200)

fig, ax = plt.subplots()
ax.plot(xin, stats.norm.pdf(xin, loc=2.5, scale=0.1))
ax.set(xlabel='beta1', ylabel='density')
None


# Compile the new Stan program.

# In[15]:


mod_prior = pystan.StanModel(model_code=stan_code_prior)


# Sample from the posterior

# In[16]:


mcmc_prior = mod_prior.sampling(data=data, chains=4, iter=2000)


# In[17]:


print(mcmc_prior)


# Plot the posterior relative to the prior.
# Note the increased error due to the prior constraints.

# In[18]:


fig, axs = plt.subplots(ncols=3, figsize=(12,6))
for ax,name in zip(axs,('beta0','beta1','sigma')):
    ax.hist(mcmc_prior[name], density=True)
    ax.axvline(true_value[name], color='black')
    ax.set(ylabel='density', title=name)
axs[1].plot(xin, stats.norm.pdf(xin, loc=2.5, scale=0.1))


# In[ ]:





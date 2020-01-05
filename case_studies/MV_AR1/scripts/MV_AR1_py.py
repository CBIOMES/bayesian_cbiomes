#!/usr/bin/env python
# coding: utf-8

# ## Multivariate autoregressive model example
# 
# In this case study we analyze multivariate time series using generic linear stochastic dynamical systems models. 
# We will use these models to analyze multispecies time series and extract information on time series interactions and forcing.
# 
# The model we will analyze is given as 
# 
# $$ \mathbf{x}(t+1) = \mathbf{\Phi} \mathbf{x}(t) + \mathbf{e}(t) $$
# 
# where $\mathbf{x}(t)$ is a multivariate time series at time $t$, $\mathbf{\Phi}$ is the linear dynamics matrix, and $\mathbf{e}(t)$ is a realization of a stochastic variable, assumed here multivariate normal with mean vector zero and covariance matrix, $\mathcal{N}(\mathbf{\mu},\mathbf{\Sigma})$.
# 
# In general, we can think of this model as a first order approximation to the dynamics of a more general multivariate dynamical system. 
# In this context, the matrix $\mathbf{\Phi}$ would represent the Jacobian matrix, containing information on the discretized partial derivatives of each state variable with respect to each of the others. 
# The matrix $\mathbf{\Sigma}$ contains information on the magnitude of stochastic forcing as well as correlations among the variables responses.

# In[1]:


import numpy as np
import pystan


# To explore this model we will first simulate some data
# 
# The procedure below uses some linear algebra (specifically the $QR$ decomposition) to generate a random matrix that will stable, i.e. that simulated time series won't blow-up.
# We do so by making the eigenvalues of the random matrix lie between 0 and 1.

# In[2]:


p = 3

SIGMA = np.diag(np.random.uniform(size=p))   #generate random numbers

eig = np.random.uniform(size=p)
Z = np.random.normal(size=(p,p))
Q,R = np.linalg.qr(Z)
d = np.diag(R)
ph = d / abs(d)
O = Q @ np.diag(ph)
PHI = O.T @ np.diag(eig) @ O


# We simulate time series given the specified parameters. 
# Notice we are drawing a multivariate normal realization every time step with mean zero and covariance matrix `SIGMA`.

# In[3]:


T = 200
y0 = np.random.uniform(size=p)
Y = np.full((p,T), fill_value=np.nan)
Y[:,0] = y0

for t in range(1,T):
    Y[:,t] = PHI @ Y[:,t-1] + np.random.multivariate_normal(mean=np.zeros(p), cov=SIGMA)


# Plot the system

# In[4]:


import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12,6))
ax.plot(Y.T)
None


# The Stan code below implements the first order multivariate model.
# As per usual, we pass Stan the data and its dimensions in the `data` block.
# In the `parameters` block, we define the linear dynamics matrix `PHI` which has as many rows and columns as there are time series.
# In this case we assume the covariance matrix is of the form
# 
# $$ \mathbf{\Sigma} = \begin{bmatrix} \sigma_{1,1} & 0 & \cdots & 0\\
#                                     0 & \sigma_{2,2} & \ddots & \vdots \\
#                                     \vdots & \ddots & \ddots & \vdots \\
#                                     0 & \cdots & \cdots & \sigma_{p,p} \end{bmatrix}   $$
#                                     
# which means that each species in the mulispecies time series responds to the environment independently.

# In[5]:


mod_code = '''data {
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
}
'''


# Compile the model

# In[6]:


mod = pystan.StanModel(model_code=mod_code)


# Organize the data for Stan
# 
# Perform mcmc

# In[7]:


data = {'p':p, 'T':T, 'Y':Y}

mcmc = mod.sampling(data=data,iter=2000,warmup=1000)

print(mcmc)


# ### Investigating the stochastic forcing matrix
# Here we relax the assumption of independent forcing and specify a covariance matrix with non zero off diagonal elements
# 
# $$ \mathbf{\Sigma} = \begin{bmatrix} \sigma_{1,1} & \sigma_{1,2} & \cdots & \sigma_{1,p}\\
#                                     \sigma_{2,1} & \sigma_{2,2} & \ddots & \vdots \\
#                                     \vdots & \ddots & \ddots & \vdots \\
#                                     \sigma_{p,1} & \cdots & \cdots & \sigma_{p,p} \end{bmatrix}   $$

# In[8]:


A = np.random.uniform(size=(p,p))*2 - 1.0
SIGMA2 = A.T @ A
print(SIGMA2)


# Simulate synthetic time series with the covariance matrix above

# In[9]:


T = 200
y20 = np.random.normal(size=p)
Y2 = np.full((p,T), fill_value=np.nan)
Y2[:,0] = y20

for t in range(1,T):
    Y2[:,t] = PHI @ Y2[:,t-1] + np.random.multivariate_normal(mean=np.zeros(p), cov=SIGMA2)


# Plot the time series

# In[10]:


fig, ax = plt.subplots(figsize=(12,6))
ax.plot(Y2.T)
None


# We estimate the full covariance in the Stan code below.
# We make use of the very handy parameter type `cov_matrix` which maintains the positive semi-definite requirement throughout the mcmc iterations.

# In[11]:


mod_code_cov = '''data {
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
}'''


# Compile the model

# In[12]:


mod_cov = pystan.StanModel(model_code=mod_code_cov)


# Package data
# 
# Perform mcmc

# In[13]:


data2 = {'p':p, 'T':T, 'Y':Y2}

mcmc_cov = mod.sampling(data=data2)

print(mcmc_cov)


# ### Incorporating informative priors on the structure of the dynamics
# 
# Here we will impose informative priors on particular elements of the $\mathbf{\Phi}$
# 
# The Stan code below imposes the following structure on $\mathbf{\Phi}$:
# 
# $$ \mathbf{\Phi} = \begin{bmatrix} \phi_{1,1} & \phi_{1,2} & 0 \\ 
#                                 \phi_{2,1} & \phi_{2,2} & \phi_{2,3}\\
#                                 0       & \phi_{3,2} & \phi_{3,3} \end{bmatrix} $$
#                                 
# The prior assumption is that $\phi_{1,3}$ and $\phi_{3,1}$ have means of zero and a very small standard deviation of $10^{-3}$

# In[14]:


mod_code_D_struc = '''data {
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
}
'''


# In[15]:


mod_struc = pystan.StanModel(model_code=mod_code_D_struc)


# In[16]:


mcmc_struc = mod_struc.sampling(data=data)
print(mcmc_struc)


import numpy as np
import pystan

#######################################################

p = 3

SIGMA = np.diag(np.random.uniform(size=p))   #generate random numbers

eig = np.random.uniform(size=p)
Z = np.random.normal(size=(p,p))
Q,R = np.linalg.qr(Z)
d = np.diag(R)
ph = d / abs(d)
O = Q @ np.diag(ph)
PHI = O.T @ np.diag(eig) @ O

#######################################################

T = 200
y0 = np.random.uniform(size=p)
Y = np.full((p,T), fill_value=np.nan)
Y[:,0] = y0

for t in range(1,T):
    Y[:,t] = PHI @ Y[:,t-1] + np.random.multivariate_normal(mean=np.zeros(p), cov=SIGMA)
	
#######################################################

import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(12,6))
ax.plot(Y.T)
None

#######################################################

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

######################################################

mod = pystan.StanModel(model_code=mod_code)

######################################################

data = {'p':p, 'T':T, 'Y':Y}

mcmc = mod.sampling(data=data,iter=2000,warmup=1000)

print(mcmc)

######################################################

A = np.random.uniform(size=(p,p))*2 - 1.0
SIGMA2 = A.T @ A
print(SIGMA2)

######################################################

T = 200
y20 = np.random.normal(size=p)
Y2 = np.full((p,T), fill_value=np.nan)
Y2[:,0] = y20

for t in range(1,T):
    Y2[:,t] = PHI @ Y2[:,t-1] + np.random.multivariate_normal(mean=np.zeros(p), cov=SIGMA2)
	
######################################################

fig, ax = plt.subplots(figsize=(12,6))
ax.plot(Y2.T)
None

######################################################

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

######################################################

mod_cov = pystan.StanModel(model_code=mod_code_cov)

######################################################

data2 = {'p':p, 'T':T, 'Y':Y2}

mcmc_cov = mod.sampling(data=data2)

print(mcmc_cov)

######################################################

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

#####################################################

mod_struc = pystan.StanModel(model_code=mod_code_D_struc)

#####################################################

mcmc_struc = mod_struc.sampling(data=data)
print(mcmc_struc)


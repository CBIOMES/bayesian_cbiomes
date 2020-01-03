import numpy as np
from scipy.integrate import odeint

#########################################################

theta = (0.25, 0.1) # γ, λ
P = 1.0
T = 30
dt = 1.0
t = np.arange(0,T+dt,dt)

##########################################################

def dxdt(P, t, γ, λ):
    return γ*P - λ*P**2
	
##########################################################

x = odeint(dxdt, P, t, args=theta)[:,0]

##########################################################

import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(t,x)
ax.set(xlabel='time', ylabel='P')
None

##########################################################

t_obs_ind = np.array((1,3,4,6,9,13,16,24))-1

obs = x[t_obs_ind] + np.random.normal(size=t_obs_ind.size, scale=0.1)
t_obs = t[t_obs_ind]

fig, ax = plt.subplots()
ax.plot(t_obs, obs, marker='o', linestyle='none')
ax.set(xlabel='time', ylabel='P observations')
None

#########################################################

stan_code = '''functions {
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
}'''

#######################################################

import pystan

mod = pystan.StanModel(model_code=stan_code)

#######################################################

data = {
    'N':len(t_obs),
    't_obs':t_obs,
    'y':obs,
    'sigma':0.1,
}

# can use chains=multiprocessing.cpu_count() to explicitly set number of chains to number of CPUs
mcmc = mod.sampling(data=data)
print(mcmc)

#######################################################

mcmc.flatnames

#######################################################

for name in 'theta','x0','x','x[1,1]':
    print('shape of {}: {}'.format(name, mcmc[name].shape))
	
#######################################################

param_names = {'theta[1]':'gamma', 'theta[2]':'lambda', 'x0':'x0'}

fig, axs = plt.subplots(ncols=3)
for ax, name in zip(axs.flat,('theta[1]','theta[2]','x0')):
    ax.hist(mcmc[name])
    ax.set(xlabel=param_names[name], ylabel='frequency')
	
#######################################################

fig, ax = plt.subplots()
ax.scatter(mcmc['theta[1]'], mcmc['theta[2]'], alpha=0.5)
ax.set(xlabel=param_names['theta[1]'], ylabel=param_names['theta[2]'])

#######################################################

theta_sin = (0.25, 0.1) 
P0 = 2.5
T_sin = 365*4
dt_sin = 1.0
t_sin = np.arange(0,T_sin+dt_sin,dt_sin)

#######################################################

def dxdt_sin(P, t, γ, λ, omega=np.pi/180.0):
    return γ*(1+np.sin(omega*t))*P - λ*P**2
	
#######################################################

x_sin = odeint(dxdt_sin, P0, t_sin, args=theta_sin)[:,0]

#######################################################

fig, ax = plt.subplots()
ax.plot(t_sin,x_sin)
ax.set(xlabel='time', ylabel='P')
None

#######################################################


t_obs_ind_sin = np.random.choice(t_sin.size, size=50, replace=False)
t_obs_ind_sin.sort()

obs_sin = np.maximum(0.0, x_sin[t_obs_ind_sin] + np.random.normal(size=t_obs_ind_sin.size, scale=0.5))
t_obs_sin = t_sin[t_obs_ind_sin]

fig, ax = plt.subplots()
ax.plot(t_obs_sin, obs_sin, marker='o', linestyle='none')
ax.set(xlabel='time', ylabel='P observations')
None

######################################################

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
}'''

###################################################

mod_sin = pystan.StanModel(model_code=stan_code_sin)

###################################################

data_sin = {
    'N':len(t_obs_sin),
    't_obs':t_obs_sin,
    'y':obs_sin,
}

mcmc_sin = mod_sin.sampling(data=data_sin, iter=2000, chains=4)
print(mcmc_sin)

####################################################

fig, axs = plt.subplots(ncols=3)
for ax, name in zip(axs.flat,('theta[1]','theta[2]','x0')):
    ax.hist(mcmc_sin[name])
    ax.set(xlabel=param_names[name], ylabel='frequency')
None

####################################################

fig, ax = plt.subplots()
ax.scatter(mcmc_sin['theta[1]'], mcmc_sin['theta[2]'], alpha=0.5)
ax.set(xlabel=param_names['theta[1]'], ylabel=param_names['theta[2]'])
None

####################################################

mu = np.mean(mcmc_sin['x'][:,:,0], axis=0)
std = np.std(mcmc_sin['x'][:,:,0], axis=0)

fig, ax = plt.subplots()
ax.fill_between(t_obs_sin, y1=mu-std, y2=mu+std, alpha=0.5)
ax.plot(t_obs_sin, mu)
ax.scatter(t_obs_sin, obs_sin)

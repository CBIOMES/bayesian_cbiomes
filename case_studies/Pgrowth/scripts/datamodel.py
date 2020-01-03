from scipy.integrate import odeint
import numpy as np


def simulate(P, t, γ, λ):
    om = 2.0*np.pi/365.0
    return γ *(1+np.sin(om*t))*P - λ*P**2
	
############################################################

def loglike(theta,obs):
    # step 1: run the model given the parameters to get state
    tobs = obs['t']
    yobs = obs['y']
    
    P = 2.5
    xhat = odeint(simulate, P, tobs, args=theta)[:,0]
    
    # step 2: compute the likelihood (compare state to observations)
    return -np.sum((yobs-xhat)**2)
	
############################################################

T = 365*4
theta = (0.25, 0.1) # γ, λ
P = 2.5
tobs = np.arange(0,T,15)
x = odeint(simulate, P, tobs, args=theta)[:,0]
yobs = x + 0.3 * np.random.normal(size=x.size)

obs = {'t':tobs, 'y':yobs}

############################################################

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(tobs, yobs)
ax.set(title='synthetic observations', xlabel='time (days)', ylabel='P')
None

############################################################

# specify the parameters
theta = (0.2,0.15) # γ, λ

# compute likelihood
L = loglike(theta,obs)
print(L)

############################################################

# gamma

gamma_grid = np.arange(0.1,0.45,0.05)
Lgamma = [loglike((g,0.1), obs) for g in gamma_grid]

fig, ax = plt.subplots()
ax.plot(gamma_grid, Lgamma, marker='o')
ax.set(title='Likelihood Profile: Gamma', xlabel='gamma', ylabel='likelihood')
None

############################################################

# lambda

lambda_grid = np.arange(0.05,0.16,0.01)
Llambda = [loglike((0.25,l), obs) for l in lambda_grid]

fig, ax = plt.subplots()
ax.plot(lambda_grid, Llambda, marker='o')
ax.set(title='Likelihood Profile: Lambda', xlabel='lambda', ylabel='likelihood')
None

############################################################

# 2D likelihood

gamma_grid = np.arange(0.2,0.31,0.01)
lambda_grid = np.arange(0.08,0.125,0.005)

Lacc = np.array([[loglike((g,l), obs) for g in gamma_grid] for l in lambda_grid])

fig, ax = plt.subplots()
cs = ax.contour(gamma_grid, lambda_grid, Lacc, levels=np.linspace(-100,0,11))
ax.clabel(cs, fontsize=10)
ax.set(title='Likelihood surface', xlabel='gamma', ylabel='lambda')
None

###########################################################

from scipy.optimize import minimize

loglike_neg = lambda theta: -loglike(tuple(theta), obs)

res = minimize(loglike_neg, (0.25,0.1), method='BFGS', options={'disp':True})
print('solution: γ={theta[0]:.3f}, λ={theta[1]:.3f}'.format(theta=res['x']))
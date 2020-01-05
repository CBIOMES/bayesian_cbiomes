#!/usr/bin/env python
# coding: utf-8

# ## Dynamics of microbial OTUs
# 
# In this notebook we use the multivariate AR(1) model to fit a multivariate time series of marine bacterial operational taxonomic units (OTUs).
# The dataset comes from
# 
# *Martin-Platero et al. 2018. High resolution time series reveals cohesive but short lived communities in coastal plankton. Nature Communications*
# 
# First load packages and read in the dataset:

# In[1]:


import pystan
import numpy as np
import pandas


# In[2]:


DAT = pandas.read_csv('data/bacterial_OTU.csv')


# Aggregate the data by phyla, which is indicated in the third column ("phylum") of the dataset:

# In[3]:


phyla = np.unique(DAT['phylum']) # extract unique phyla IDs

# sum all OTUs of that phylum
# use sums as rows for new DataFrame
PHY = pandas.DataFrame([DAT.loc[DAT['phylum']==phylum].iloc[:,8:].sum() for phylum in phyla])
# associate each row with phylum
PHY = PHY.set_index(phyla)


# Now take a look to see how each phyla contributes to the total abundances:

# In[4]:


PHY.sum(axis=1)


# In[5]:


index_max4 = np.argsort(PHY.sum(axis=1).values)[-4:][::-1]
# for consistency with the R code, rearrange the index (now they are in order of appearance in dataset)
index_max4 = index_max4[np.array((1,0,3,2))]

PHY = PHY.iloc[index_max4]
PHY


# In[6]:


import matplotlib.pyplot as plt

fig,ax = plt.subplots(figsize=(12,6))
ax.plot(PHY.values.T)
ax.legend(PHY.index)
None


# In[7]:


dat_PHY = {'T':PHY.shape[1], 'p':PHY.shape[0], 'Y':PHY.values}


# In[8]:


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
    Y[,1] ~ normal(init, sigma);            //distribution of the initial conditions
    for(i in 2:T){
        Y[,i] ~ normal(PHI*Y[,i-1],sigma);  //conditional predictive distribution
    }
}
'''


# In[9]:


mod = pystan.StanModel(model_code=mod_code)


# In[10]:


mcmc = mod.sampling(data=dat_PHY, iter=2000, warmup=1000)
print(mcmc)


# In[ ]:





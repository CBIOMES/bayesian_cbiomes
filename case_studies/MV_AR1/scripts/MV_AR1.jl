using StanSample, Statistics, PyPlot, Random, Distributions, LinearAlgebra

###################################################

p = 3;
SIGMA = diagm(rand(Uniform(),p))

eig = rand(Uniform(),p)
Z = rand(Normal(),p,p)
Q,R = qr(Z)
d = diag(R)
ph = d ./ abs.(d)
O = Q * diagm(ph)
PHI = transpose(O) * diagm(eig) * O;

###################################################

T = 200
y0 = rand(Normal(),p)
Y = zeros(p,T)
Y[:,1] = y0

for i in 2:T
    Y[:,i] = PHI * Y[:,i-1] + rand(MvNormal(zeros(3),SIGMA))
end

###################################################

fig, ax = plt.subplots(figsize=(12,6))
ax.plot(transpose(Y));

###################################################

mod_code = "data {
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
}";

#################################################

data = Dict("p" => p, "T" => T, "Y" => Y)

#################################################

sm = SampleModel("MV_AR", mod_code)

#################################################

(sample_file, log_file) = stan_sample(sm, data=data, n_chains = 4);

#################################################

chns = read_samples(sm)

#################################################

ESS = ess(chns)

#################################################

A = rand(Uniform(),p,p) .* 2 .- 1
SIGMA2 = transpose(A) * A

#################################################

T = 200
y20 = rand(Normal(),p)
Y2 = zeros(p,T)
Y2[:,1] = y20

for i in 2:T
    Y2[:,i] = PHI * Y2[:,i-1] + rand(MvNormal(zeros(3),SIGMA2))
end

#################################################

fig, ax = plt.subplots(figsize=(12,6))
ax.plot(transpose(Y2));

#################################################


mod_code_cov = "data {
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
}";

################################################

data2 = Dict("p" => p, "T" => T, "Y" => Y2)

################################################

sm2 = SampleModel("MV_AR2", mod_code_cov)

################################################

(sample_file, log_file) = stan_sample(sm2, data=data2, n_chains = 4);

################################################

chns2 = read_samples(sm2)

################################################

ESS2 = ess(chns2)

################################################

mod_code_D_struc = "data {
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
}";

#################################################

sm3 = SampleModel("MV_AR3", mod_code_D_struc)

#################################################

(sample_file, log_file) = stan_sample(sm3, data=data, n_chains = 4);

#################################################

chns3 = read_samples(sm3)

#################################################

ESS3 = ess(chns3)



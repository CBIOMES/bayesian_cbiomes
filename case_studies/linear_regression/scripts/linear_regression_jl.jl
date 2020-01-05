using StanSample, CSV, DataFrames, PyPlot, Statistics, Distributions, Random

versioninfo()

N = 100;
beta1 = 1.5;
beta0 = 2.0;
sigma = 1.25;

x = rand(Normal(),N);
y = beta0 .+ beta1 .* x .+ rand(Normal(0,sigma),N);

fig, ax = PyPlot.subplots()
ax.scatter(x,y)
ax.set_ylabel("y",fontsize=16)
ax.set_xlabel("x",fontsize=16);

stancode = "
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
}";

mod = SampleModel("linear_reg", stancode)

data = Dict("N" => N, "x" => x, "y" => y);

(sample_file, log_file) = stan_sample(mod, data=data, n_chains = 4);

chns = read_samples(mod)

ESS = ess(chns)

rawdata = DataFrame(chns, showall=true, sorted=true, append_chains=true);

cnames=["beta0","beta1","sigma"]
fig, axs = PyPlot.subplots(3, 1, figsize = (5,6))
for i in 1:3
    axs[i].scatter(collect(1:1:4000), rawdata[:,i],s=1)
    axs[i].plot(collect(1:1:4000),zeros(4000).+chns.info[1].x[2][1][i,2],"r")
    axs[i].set_title(cnames[i], ha="center", fontsize=12, color = "k");
end
fig.subplots_adjust(bottom=0.1, top=0.96, left=0.1, right=0.95,
                    wspace=0.2, hspace=0.4)

fig, axs = PyPlot.subplots(1, 3, figsize = (8,2))
for i in 1:3
    axs[i].hist(rawdata[:,i],bins=20)
    axs[i].axvline(chns.info[1].x[2][1][i,2],color="r")
    axs[i].set_title(cnames[i], ha="center", fontsize=12, color = "k");
end

fig.subplots_adjust(bottom=0.1, top=0.96, left=0.1, right=0.95,
                    wspace=0.4, hspace=0.1)

beta0_mean = chns.info[1].x[2][1][1,2];
beta1_mean = chns.info[1].x[2][1][2,2];
sigma_mean = chns.info[1].x[2][1][3,2];
fig, ax = PyPlot.subplots()
ax.scatter(x,y)
for i in 1:20
    ax.plot(x,rawdata[i,1] .+ x .* rawdata[i,2],"k--",alpha=0.2,lw=0.5)
end
ax.plot(x,beta0_mean .+ x .* beta1_mean,"r")
ax.set_ylabel("y",fontsize=16)
ax.set_xlabel("x",fontsize=16);

stancode_prior = "
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
}";

xin = collect(0:0.01:3);

using StatsPlots

StatsPlots.plot(xin,Normal(2.5,0.1),xlabel = "beta1",ylabel = "Density")

mod_prior = SampleModel("LRprior", stancode_prior)

(sample_file, log_file) = stan_sample(mod_prior, data=data, n_chains = 4);

chns_prior = read_samples(mod_prior)

ESS_prior = ess(chns_prior)

rawdata_prior = DataFrame(chns_prior, showall=true, sorted=true, append_chains=true);

fig, axs = PyPlot.subplots(1, 3, figsize = (8,2))
for i in 1:3
    axs[i].hist(rawdata_prior[:,i],bins=20)
    axs[i].axvline(chns_prior.info[1].x[2][1][i,2],color="r")
    axs[i].set_title(cnames[i], ha="center", fontsize=12, color = "k");
end
axs[1].axvline(beta0,color="m")
axs[2].axvline(beta1,color="m")
axs[3].axvline(sigma,color="m")
axs[2].set_xlim((0.4,3))
fig.subplots_adjust(bottom=0.1, top=0.96, left=0.1, right=0.95,
                    wspace=0.4, hspace=0.1)



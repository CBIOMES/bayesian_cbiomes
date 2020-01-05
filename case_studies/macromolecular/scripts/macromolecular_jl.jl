using StanSample, CSV, DataFrames, PyPlot, Statistics

versioninfo()

data = CSV.read("data/macromolecules.csv");
NH4 = CSV.read("data/flynn_macromolecules.csv");

dat = join(data,NH4,on=:t);
dat = hcat(dat[:,1:4],dat[:,7]);
rename!(dat,Dict(:x1 => "NH4"));
y_data = zeros(size(dat,1),4)
y_data[:,1] = dat[:,4] .* 0.001;
y_data[:,2] = dat[:,2] .* 0.001;
y_data[:,3] = dat[:,3] .* 1000;
y_data[:,4] = dat[:,5] .* 0.08325909;

fig, axs = PyPlot.subplots(2, 2, figsize = (6,4))
axs[1].scatter(dat[:,1], y_data[:,1])
axs[2].scatter(dat[:,1], y_data[:,2])
axs[3].scatter(dat[:,1], y_data[:,3])
axs[4].scatter(dat[:,1], y_data[:,4])
axs[1].set_title("C Quota")
axs[2].set_title("N Quota")
axs[3].set_title("Chl Quota")
axs[4].set_title("Ammonium")
fig.subplots_adjust(wspace=0.3, hspace=0.44)

###### macromolecular model stan code
const macromolecularmodel = "functions {
  real[] macro(real   t,           // time
               real[] x,           // state x[1]:CH  x[2]:PR, x[3]:Chl , x[4]:N
               real[] theta,
               real[] x_r,
               int[]  x_i) {       // parameters

    real CNpro = theta[1]; 
    real KN    = theta[2];    
    real mu    = theta[3]; 
    real CHsyn = theta[4]; 
    real m_ex  = theta[5];  
    real R_ex  = theta[6];  
    real tau   = theta[7];
    real b     = theta[8];

    real PRsynth = theta[3]*x[4]/(theta[2]+x[4]);
    real r0      = theta[8]*(x[2]/x[1]);
    real Chl     = x[3]*x[2];
    real Rcell   = x[1]/x[2];
    real excr    = (1/2)*theta[5]*(1+tanh(Rcell - theta[6]));
    
    real dCH    = x[2]*(theta[4] - excr);
    real dr     = (1/theta[7])*(r0-x[3]);
    real dPR    = x[2]*PRsynth;
    real dN     = -dPR/(1+exp(-10000*x[4])); 

    return {dCH,dPR,dr,dN};
  }
}
data {
  int<lower = 0> n;           // num obs
  real t_obs[n];              // obs times
  real<lower = 0> y[n,4];     // observed variable at measurement times
}
parameters {
  real<lower = 0> theta[8];   // parameters
  real<lower = 0> x0[4];      // initial population
  real<lower = 1E-15> sigma[4]; // obs error
}
transformed parameters {
  real x[n,4] = integrate_ode_rk45(macro, x0, -1, t_obs, theta, rep_array(0.0,0), rep_array(0,0), 1e-6, 1e-5, 1e3) ;
  for(i in 1:n){
    x[i,3] = x[i,3]*x[i,2]*1E6;
  }
}
model {
  x0[1]    ~ normal(0.1,1);
  x0[2]    ~ normal(0.1,1);
  x0[3]    ~ normal(10,10);
  x0[4]    ~ normal(0.1,1);  
  theta[1] ~ normal(6.6,10); //prior on CHpro
  theta[2] ~ normal(0.002,3); //prior on KN
  theta[3] ~ normal(0.3,1);
  theta[4] ~ normal(5,10);
  theta[5] ~ normal(10,10);
  theta[6] ~ normal(13,10);
  theta[7] ~ normal(10,10); 
  theta[8] ~ normal(0.05,5); 
  
  for(i in 1:4){
    y[1:n,i] ~ normal(x[1:n,i], sigma[i]);
  }
}";

macromoleculardata = Dict("n" => size(dat,1), "t_obs" => dat[:,1], "y" => y_data);

sm = SampleModel("MacroMolecularModel", macromolecularmodel)

(sample_file, log_file) = stan_sample(sm, data=macromoleculardata, n_chains = 4);

chns = read_samples(sm)

ESS = ess(chns)

rawdata = DataFrame(chns, showall=true, sorted=true, append_chains=true);

rawdata[:,17:end]
dat_summary = zeros(35,4,4)
for i in 1:35
    for j in 1:4
        dat_summary[i,j,1] = mean(rawdata[:,16+4*(i-1)+j])
        dat_summary[i,j,2] = std(rawdata[:,16+4*(i-1)+j])
        dat_summary[i,j,3] = quantile(rawdata[:,16+4*(i-1)+j],0.025)
        dat_summary[i,j,4] = quantile(rawdata[:,16+4*(i-1)+j],0.975)
    end
end

dat_mean = chns.info[1].x[2][1][17:end,2];
dat_std = chns.info[1].x[2][1][17:end,3];
dat_97 = chns.info[1].x[2][2][17:end,6];
dat_2 = chns.info[1].x[2][2][17:end,2];

data_mean = reverse(rotl90(reshape(dat_mean,4,35)),dims=1)
data_std = reverse(rotl90(reshape(dat_std,4,35)),dims=1)
data_97 = reverse(rotl90(reshape(dat_97,4,35)),dims=1)
data_2 = reverse(rotl90(reshape(dat_2,4,35)),dims=1);

fig, axs = PyPlot.subplots(2, 2, figsize = (6,4))
axs[1].scatter(dat[:,1], y_data[:,1])
axs[1].plot(dat[:,1], data_mean[:,1],"r")
axs[1].plot(dat[:,1], data_2[:,1],"r--")
axs[1].plot(dat[:,1], data_97[:,1],"r--")
axs[2].scatter(dat[:,1], y_data[:,2])
axs[2].plot(dat[:,1], data_mean[:,2],"r")
axs[2].plot(dat[:,1], data_2[:,2],"r--")
axs[2].plot(dat[:,1], data_97[:,2],"r--")
axs[3].scatter(dat[:,1], y_data[:,3])
axs[3].plot(dat[:,1], data_mean[:,3],"r")
axs[3].plot(dat[:,1], data_2[:,3],"r--")
axs[3].plot(dat[:,1], data_97[:,3],"r--")
axs[4].scatter(dat[:,1], y_data[:,4])
axs[4].plot(dat[:,1], data_mean[:,4],"r")
axs[4].plot(dat[:,1], data_2[:,4],"r--")
axs[4].plot(dat[:,1], data_97[:,4],"r--")
axs[1].set_title("C Quota")
axs[2].set_title("N Quota")
axs[3].set_title("Chl Quota")
axs[4].set_title("Ammonium")
fig.subplots_adjust(wspace=0.3, hspace=0.44)

fig, axs = PyPlot.subplots(2, 2, figsize = (6,4))
axs[1].scatter(dat[:,1], y_data[:,1])
axs[1].plot(dat[:,1], dat_summary[:,1,1],"r")
axs[1].plot(dat[:,1], dat_summary[:,1,3],"r--")
axs[1].plot(dat[:,1], dat_summary[:,1,4],"r--")
axs[2].scatter(dat[:,1], y_data[:,2])
axs[2].plot(dat[:,1], dat_summary[:,2,1],"r")
axs[2].plot(dat[:,1], dat_summary[:,2,3],"r--")
axs[2].plot(dat[:,1], dat_summary[:,2,4],"r--")
axs[3].scatter(dat[:,1], y_data[:,3])
axs[3].plot(dat[:,1], dat_summary[:,3,1],"r")
axs[3].plot(dat[:,1], dat_summary[:,3,3],"r--")
axs[3].plot(dat[:,1], dat_summary[:,3,4],"r--")
axs[4].scatter(dat[:,1], y_data[:,4])
axs[4].plot(dat[:,1], dat_summary[:,4,1],"r")
axs[4].plot(dat[:,1], dat_summary[:,4,3],"r--")
axs[4].plot(dat[:,1], dat_summary[:,4,4],"r--")
axs[1].set_title("C Quota")
axs[2].set_title("N Quota")
axs[3].set_title("Chl Quota")
axs[4].set_title("Ammonium")
fig.subplots_adjust(wspace=0.3, hspace=0.44)

cnames=["CNpro","KN","mu","CHsyn","m_ex","R_ex","tau","b"]
fig, axs = PyPlot.subplots(4, 2, figsize = (8,8))
for i in 1:8
    axs[i].scatter(collect(1:1:4000), rawdata[:,4+i],s=1)
    axs[i].plot(collect(1:1:4000),zeros(4000).+chns.info[1].x[2][1][4+i,2],"r")
    axs[i].set_title(cnames[i], ha="center", fontsize=12, color = "k");
end
fig.subplots_adjust(bottom=0.1, top=0.96, left=0.1, right=0.95,
                    wspace=0.2, hspace=0.4)

fig, axs = PyPlot.subplots(4, 2, figsize = (8,8))
for i in 1:8
    axs[i].hist(rawdata[:,4+i],bins=20)
    axs[i].axvline(chns.info[1].x[2][1][4+i,2],color="r")
    axs[i].set_title(cnames[i], ha="center", fontsize=12, color = "k");
end

fig.subplots_adjust(bottom=0.1, top=0.96, left=0.1, right=0.95,
                    wspace=0.2, hspace=0.4)



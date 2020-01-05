using DifferentialEquations, StanSample, PyPlot, DataFrames, Statistics

# specify indices
i_n = 1
i_p = 2
i_z = 3

# specify parameter values
#        vmax,  nuthalfsat, graz, mort_p, mort_z, irr
theta = [0.075, 0.3,        0.02, 0.02,   0.03,   0.8]

# vmax:       maximum growth rate in Michaelis Menten formulation
# nuthalfsat: nutrient half saturation in Michaelis Menten formulation
# graz:       zooplankton grazing rate
# mort_p:     phytoplankton mortality rate
# mort_z:     zooplankton mortality rate
# irr:        light amplitude


# specify model 
function dxdt!(du,u,p,t)   
    light = 1.0 + 0.5*(p[6]*sin(π*((t-81.25)/182.5))-p[6])
    growth = p[1]*u[i_n]/(p[2]+u[i_n]) * light * u[i_p]
    grazing = p[3]*u[i_p]*u[i_z]
    ploss = p[4]*u[i_p]
    zloss = p[5]*u[i_z]*u[i_z]
    du[1] = -growth+ploss+zloss # N
    du[2] = growth-grazing-ploss # P
    du[3] = grazing-zloss       # Z
end

# initial conditions
x0 = [0.6,0.15,0.23]

# initialize time vector
t = collect(0.0:1.0:2*365.0);

tspan = (0.0,2*365.0)

prob = ODEProblem(dxdt!, x0, tspan, theta)
sol = solve(prob,reltol=1e-6,saveat=1.0);

# reorganize sol.u into 2d array
x = reverse(rotl90(hcat(sol.u...)),dims=1);

iobs = sample(1:size(t,1)÷10,20,replace=false)*10;
iobs = sort(iobs)
tobs = t[iobs]
iobsvar = [i_p,i_z]
sigma = [0.03 0.03]
obs = max.(0.0, x[iobs,iobsvar] .+ randn(20,2) .* sigma);

colors = ("#1f77b4","#2ca02c","#d62728") # colors for N, P, Z

light = 1.0 .+ 0.5.*(theta[6] .* sin.(π.*((t.-81.25)./182.5)) .- theta[6])
cname = ["N", "P", "Z"]
fig,ax = plt.subplots(figsize=(12,10))
for i in 1:3
    ax.plot(t,x[:,i], color=colors[i], label=cname[i])
end
    ax.plot(t,theta[1]*light,linestyle=":",color="#ff7f0e",label="light (scaled)")
for i in 1:2
    ax.plot(tobs, obs[:,i], color=colors[i+1], marker="o", ls="none")
end
ax.set(title="NPZ reference simulation and synthetic observations", xlabel="time (days)")
ax.legend();

stan_code = """
functions {
   real[] npz(real t,       // time
              real[] x,     // state
              real[] theta, // parameters
              real[] x_r,   // fixed real data (empty)
              int[] x_i) {  // fixed integer data (empty)
   
    /*
    guide to theta:
    theta[1]:  vmax         maximum growth rate in Michaelis Menten formulation
    theta[2]:  nuthalfsat   nutrient half saturation in Michaelis Menten formulation
    theta[3]:  graz         zooplankton grazing rate
    theta[4]:  mort_p       phytoplankton mortality rate
    theta[5]:  mort_z       zooplankton mortality rate
    theta[6]:  irr          light amplitude
    */

    real n = fmax(0.0, x[1]);
    real p = fmax(0.0, x[2]);
    real z = fmax(0.0, x[3]);

    real light = 1.0 + 0.5*(theta[6]*sin(pi()*((t-81.25)/182.5))-theta[6]); 
    real growth = theta[1]*n/(theta[2]+n) * light * p;
    real grazing = theta[3]*p*z;
    real ploss = theta[4]*p;
    real zloss = theta[5]*z*z;
    
    return {-growth+ploss+zloss,growth-grazing-ploss,grazing-zloss};
  }
}
data {
    int<lower=0> nobs;               // number of timesteps with observations
    real tobs[nobs];                 // obs times
    int<lower=0> nobsvar;            // number of observed variables
    int<lower=0> iobsvar[nobsvar];   // index of observed variable (N=1, P=2, Z=3)
    real<lower=0> obs[nobs,nobsvar]; // observed variable at measurement times
}
parameters {
    real<lower=0> vmax;
    real<lower=0> nuthalfsat;
    real<lower=0> graz;
    real<lower=0> mort_p;
    real<lower=0> mort_z;
    real<lower=0,upper=1> irr;
    real<lower=0> x0[3];            // initial conditions
    real<lower=0> sigma[nobsvar];   // obs error
}
transformed parameters {
    real theta[6] = {vmax,nuthalfsat,graz,mort_p,mort_z,irr};
    real x[nobs, 3] = integrate_ode_rk45(npz, x0, 0, tobs, theta,
                                         rep_array(0.0, 0), rep_array(0, 0),
                                         1e-5, 1e-4, 1e4);
}
model {
    vmax       ~ normal(0.1, 0.1);
    nuthalfsat ~ uniform(0.0, 1.0);
    graz       ~ normal(0.01, 0.01);
    mort_p     ~ normal(0.01, 0.01);
    mort_z     ~ normal(0.01, 0.01);
    irr        ~ uniform(0.0, 1.0);
    x0[1:3]    ~ normal(0.1, 0.1);
    for (iobs in 1:nobs){
        obs[iobs,] ~ normal(x[iobs,iobsvar], sigma);
    }
}
""";

sm = SampleModel("NPZ", stan_code)

data = Dict("nobs" => size(tobs,1), "tobs" => tobs, "nobsvar" => size(iobsvar,1), "iobsvar" => iobsvar, "obs" => obs)

(sample_file, log_file) = stan_sample(sm, data=data, n_chains = 4);

chns = read_samples(sm)

ESS = ess(chns)

rawdata = DataFrame(chns, showall=true, sorted=true, append_chains=true);

# reshape the results
dat_summary = zeros(20,3,4)
for i in 1:20
    for j in 1:3
        dat_summary[i,j,1] = mean(rawdata[:,17+3*(i-1)+j])
        dat_summary[i,j,2] = quantile(rawdata[:,17+3*(i-1)+j],0.25)
        dat_summary[i,j,3] = quantile(rawdata[:,17+3*(i-1)+j],0.5)
        dat_summary[i,j,4] = quantile(rawdata[:,17+3*(i-1)+j],0.75)
    end
end

fig,ax = plt.subplots(figsize=(12,10))

for i in 1:3
    ax.plot(t,x[:,i], color=colors[i], label=cname[i])
    ax.errorbar(x=tobs, y=dat_summary[:,i,3], yerr=[dat_summary[:,i,3]-dat_summary[:,i,2],dat_summary[:,i,4]-dat_summary[:,i,3]], 
        ls=":", color=colors[i], marker="+", capsize=1.0)
end
for i in 1:2
    ax.plot(tobs, obs[:,i], color=colors[i+1], marker="o", ls="none")
end
ax.set(title="NPZ reference simulation, synthetic observations and model posterior at obs locations", xlabel="time (days)")
ax.legend();

# reshape rawdata
data_bp = zeros(4000,20,3)
for i in 1:20
    for j in 1:3
        data_bp[:,i,j] = rawdata[:,17+3*(i-1)+j]
    end
end

positions = collect(1:1:20)
fig,axs = plt.subplots(nrows=3, sharex=true, figsize=(12,10))
for i in 1:3
    axs[i].boxplot(data_bp[:,:,i], positions=positions, widths=0.7, medianprops=Dict("color"=>colors[i]),
    flierprops=Dict("marker"=>".", "alpha"=>0.5))
    axs[i].plot(positions, x[iobs,i], marker="x", markeredgewidth=3, markersize=10, ls="none", label="true state", color=colors[i])
    axs[i].set(ylabel=cname[i])
end
for i in 1:2
    axs[i+1].plot(positions, obs[:,i], marker="o", markersize=10, ls="none", label="observation", color=colors[i+1])
end
for ax in axs
    ax.legend()
end
axs[1].set(title="Posterior model state at the observation locations")
axs[3].set(xticklabels=tobs, xlabel="time (days)");

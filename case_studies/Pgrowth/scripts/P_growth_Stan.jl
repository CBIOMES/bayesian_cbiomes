using DifferentialEquations, PyPlot, Distributions, LaTeXStrings, StanSample

#########################################################################

theta = [0.25, 0.1] # γ, λ
P = 1.0
T = 30.0
dt = 1.0
t = collect(0:dt:T);

##########################################################################

dxdt(u,p,t) = p[1]*u - p[2]*u^2
tspan = (0,T)
prob = ODEProblem(dxdt, P, tspan, theta)

##########################################################################

sol = solve(prob,reltol=1e-6,saveat=1.0);

##########################################################################

fig, ax = PyPlot.subplots()
ax.plot(sol.u)
ax.set_ylabel("P",fontsize=16)
ax.set_xlabel("time",fontsize=16);

##########################################################################

t_obs_ind = [1,3,4,6,9,13,16,24];
obs = sol.u[t_obs_ind] .+ rand(Normal(0,0.1),size(t_obs_ind,1))
t_obs = sol.t[t_obs_ind];
fig, ax = PyPlot.subplots()
ax.scatter(t_obs,obs)
ax.set_ylabel("P observation",fontsize=16)
ax.set_xlabel("time",fontsize=16);

##########################################################################

stancode = "functions {
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
}";

##########################################################################

pgrowthdata = Dict("N" => size(t_obs,1), "t_obs" => t_obs, "y" => obs, "sigma" => 0.1);

##########################################################################

sm = SampleModel("PgrowthModel", stancode)

##########################################################################

(sample_file, log_file) = stan_sample(sm, data=pgrowthdata, n_chains = 4);

##########################################################################

chns = read_samples(sm)

##########################################################################

ESS = ess(chns)

##########################################################################

rawdata = DataFrame(chns, showall=true, sorted=true, append_chains=true);

##########################################################################

cnames=[L"\gamma",L"\lambda","x0"]
fig, axs = PyPlot.subplots(1, 3, figsize = (8,2))
for i in 1:3
    axs[i].hist(rawdata[:,i],bins=20)
    axs[i].axvline(chns.info[1].x[2][1][i,2],color="r")
    axs[i].set_title(cnames[i], ha="center", fontsize=12, color = "k");
end

fig.subplots_adjust(bottom=0.1, top=0.96, left=0.1, right=0.95,
                    wspace=0.4, hspace=0.1)
					
#########################################################################

fig, ax = PyPlot.subplots()
ax.scatter(rawdata[:,1],rawdata[:,2],marker = "o")
ax.set_ylabel(L"\lambda",fontsize=16)
ax.set_xlabel(L"\gamma",fontsize=16);

#########################################################################

theta_sin = [0.25, 0.1]
P0 = 2.5
T_sin = 365*4
dt_sin = 1.0
t_sin = collect(0:dt_sin:T_sin);

#########################################################################

dxdt_sin(u,p,t) = p[1]*(1+sin(2*π/365*t))*u - p[2]*u^2
tspan_sin = (0.0,T_sin)
prob_sin = ODEProblem(dxdt_sin, P0, tspan_sin, theta_sin)

#########################################################################

sol_sin= solve(prob_sin,reltol=1e-6,saveat=1.0);

#########################################################################

fig, ax = PyPlot.subplots()
ax.plot(sol_sin.u)
ax.set_ylabel("P",fontsize=16)
ax.set_xlabel("time",fontsize=16);

#########################################################################

t_obs_ind_sin = sample(1:size(t_sin,1),50,replace=false);
t_obs_ind_sin = sort(t_obs_ind_sin);
t_obs_sin = t_sin[t_obs_ind_sin]
obs_sin = sol_sin.u[t_obs_ind_sin].+ rand(Normal(0,0.1),size(t_obs_ind_sin,1));

#########################################################################

fig, ax = PyPlot.subplots()
ax.scatter(t_obs_sin,obs_sin)
ax.set_ylabel("P observation",fontsize=16)
ax.set_xlabel("time",fontsize=16);

#########################################################################

stancode_sin = "functions {
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
}";

#########################################################################

pgrowthdata_sin = Dict("N" => size(t_obs_sin,1), "t_obs" => t_obs_sin, "y" => obs_sin);

#########################################################################

sm_sin = SampleModel("PgrowthModel_sin", stancode_sin)

#########################################################################

(sample_file, log_file) = stan_sample(sm_sin, data=pgrowthdata_sin, n_chains = 4);

#########################################################################

chns_sin = read_samples(sm_sin)

#########################################################################

ESS = ess(chns_sin)

#########################################################################

rawdata_sin = DataFrame(chns_sin, showall=true, sorted=true, append_chains=true);

#########################################################################


cnames=[L"\gamma",L"\lambda","x0"]
fig, axs = PyPlot.subplots(1, 3, figsize = (8,2))
for i in 1:3
    axs[i].hist(rawdata_sin[:,i+1],bins=20)
    axs[i].axvline(chns_sin.info[1].x[2][1][i+1,2],color="r")
    axs[i].set_title(cnames[i], ha="center", fontsize=12, color = "k");
end

fig.subplots_adjust(bottom=0.1, top=0.96, left=0.1, right=0.95,
                    wspace=0.4, hspace=0.1)
					
########################################################################

fig, ax = PyPlot.subplots()
ax.scatter(rawdata_sin[:,2],rawdata_sin[:,3],marker = "o")
ax.set_ylabel(L"\lambda",fontsize=16)
ax.set_xlabel(L"\gamma",fontsize=16);

########################################################################

fig, ax = PyPlot.subplots()
ax.scatter(t_obs_sin,obs_sin,color="r")
ax.plot(t_obs_sin,chns_sin.info[1].x[2][1][5:end,2])
ax.set_ylabel("P observation",fontsize=16)
ax.set_xlabel("time",fontsize=16);

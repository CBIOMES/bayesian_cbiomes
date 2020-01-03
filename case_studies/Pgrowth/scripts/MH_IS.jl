using DifferentialEquations, PyPlot, Distributions, LaTeXStrings, Statistics, CSV

#####################################################

theta = [0.25, 0.1] # γ, λ
T = 365*4
t = collect(0:1:365*4)
P0 = 1.0;

#####################################################

simulate(u,p,t) = p[1]*(1+sin(2*π/365*t))*u - p[2]*u^2

#####################################################

function loglike(theta, obs)
    u₀ = 2.5;
    tspan = (0.0, T)
    p = [theta[1],theta[2]]
    prob = ODEProblem(simulate,u₀,tspan,p)
    xhat = solve(prob,reltol=1e-6,saveat=15.0).u
    L = -sum((obs[:,2] .- xhat).^2)
    return L
end

######################################################

data = CSV.read("data/obs.csv")
obs = hcat(data[:,1],data[:,2]);

######################################################

# setup
nc = 1000 # number of iterations (length of the chain)
theta_acc = zeros(nc,2) # accumulator for sample (parameter posterior)

thetaold = [0.2, 0.1] # initial condition for chain (parameter values)

# priors (use uniforms so set upper and lower bounds)
gam_low = 0.2
gam_up = 0.3
lam_low = 0.09
lam_up = 0.11;

######################################################

for i in 1:nc
    # Step 1: Generate Trial Candidate (draw from parameter prior)
    gammac = rand(Uniform(gam_low, gam_up))
    lambdac = rand(Uniform(lam_low, lam_up))
    thetac = [gammac,lambdac]

    # Step 2: Compute acceptance probability for this candidate
    # compute log likelihood for candidate
    like_num = loglike(thetac,obs)

    # compute log likelihood for last member of chain (last sample member)
    like_denom = loglike(thetaold,obs)

    # compute the likelihood ratio
    likeratio = exp(like_num - like_denom)

    # determine acceptance probability
    A = min(1.0,likeratio) 

    # Step 3: Choose whether or not to accept candidate
    # accept candidate with probability A, otherwise revert to previous chain member
    if rand(Uniform()) ≤ A
        thetanew = thetac
    else
        thetanew = thetaold
    end
    
    # add to the sample 
    theta_acc[i,:] = thetanew
    thetaold = thetanew
end

######################################################

fig, axs = plt.subplots(nrows=2, sharex=true)
axs[1].plot(theta_acc[:,1])
axs[1].axhspan(gam_low, gam_up, color="green", alpha=0.3)
axs[1].set(ylabel=L"\gamma", title="trace plots")
axs[2].plot(theta_acc[:,2])
axs[2].axhspan(lam_low, lam_up, color="green", alpha=0.3)
axs[2].set(ylabel=L"\lambda", xlabel="chain iteration");

######################################################

using KernelDensity

######################################################

cname = [L"\gamma", L"\lambda"]
fig, axs = plt.subplots(nrows=2, ncols=2, sharex="col")
for icol in 1:2
    # histogram
    axs[1,icol].hist(theta_acc[:,icol])
    # kernel smooth density
    ksd = kde(theta_acc[:,icol])
    x = collect(minimum(theta_acc[:,icol]):0.0001:maximum(theta_acc[:,icol]))
    axs[2,icol].plot(x,pdf(ksd,x))
    axs[2,icol].set_xlabel(cname[icol])
    if icol == 1
        axs[1,icol].set_ylabel("frequency")
        axs[2,icol].set_ylabel("density")
    end
end

#####################################################

# scatter plot of joint distribution
fig, ax = plt.subplots()
ax.scatter(theta_acc[:,1], theta_acc[:,2])
ax.set(title="scatter plot", xlabel=L"\gamma", ylabel=L"\lambda");

#####################################################


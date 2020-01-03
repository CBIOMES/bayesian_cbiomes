using StanSample, Statistics, PyPlot, CSV, DataFrames, LinearAlgebra

##############################################################

DAT = CSV.read("data/bacterial_OTU.csv");

##############################################################

phyla = unique(DAT[:,3])
PHY = zeros(length(phyla),size(DAT[:,9:end],2))
for i in 1:length(phyla)
    PHY[i,:] = Matrix(aggregate(DAT[findall(x -> x == phyla[i], DAT.phylum),9:end],sum))
end
PHY = hcat(phyla,PHY)
PHY_sum = convert(DataFrame,hcat(phyla,sum(PHY[:,2:end],dims=2)))
sort!(PHY_sum,:x2,rev=true)

##############################################################

phy = zeros(4,size(PHY,2)-1)
for i in 1:4
    phy[i,:] = PHY[findall(x->x==PHY_sum[i,1],PHY)[1][1],2:end]
end

##############################################################

fig,ax = plt.subplots(figsize=(12,6))
ax.plot(transpose(phy))
ax.legend(PHY_sum[1:4,1]);

###############################################################

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
    Y[,1] ~ normal(init, sigma);            //distribution of the initial conditions
    for(i in 2:T){
        Y[,i] ~ normal(PHI*Y[,i-1],sigma);  //conditional predictive distribution
    }
}";

##############################################################

dat_PHY = Dict("T" => size(phy,2),"p" => size(phy,1), "Y" => phy)

##############################################################

sm = SampleModel("MV_AR_OTU", mod_code)

##############################################################

(sample_file, log_file) = stan_sample(sm, data=dat_PHY, n_chains = 4);

##############################################################

chns = read_samples(sm)

##############################################################

ESS = ess(chns)




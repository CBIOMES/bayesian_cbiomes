library(rstan)
options(mc.cores=parallel::detectCores())

DAT <- read.csv('data/bacterial_OTU.csv',stringsAsFactors=FALSE)

phyla <- unique(DAT[,3])   #extract unique phyla IDS
PHY   <- data.frame()      #open empty data frame
for(i in 1:length(phyla)){
	xtmp <- apply(as.data.frame(DAT[DAT[,3]==phyla[i],9:ncol(DAT)]),2,sum) #sum all OTUs of that phyla
	PHY  <- rbind(PHY,xtmp)                                                #attach as rows to the empty data frame
}

rbind(1:nrow(PHY),rowSums(PHY))  #list row number alongside row sums

phy <- PHY[1:4,]

options(repr.plot.width=6, repr.plot.height=4)
matplot(t(phy),type='l')

dat_PHY <- list(T=ncol(phy),
                p=nrow(phy),
                Y=phy)

mod_code <- "data {
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
}"

mod <- stan_model(model_code=mod_code)

mcmc <- sampling(mod,data=dat_PHY,iter=2000,warmup=1000,open_progress=TRUE)

mcmc

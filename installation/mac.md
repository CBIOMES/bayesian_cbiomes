# Instructions to install Jupyter notebooks
*NOTE: If you encounter any difficulty with these instructions, please create an issue in GitHub. 
We can help you through the installation if you get stuck and would also like to hear about issues even if you fixed them yourself.*

The first step is to install Anaconda via the instructions here: https://docs.anaconda.com/anaconda/install/.
You can download Anaconda 2 or Anaconda 3.
You will get Jupyter notebooks along with that download.


# Instructions to install Stan
Instructions to install Stan on Max will depend on whether you are using R (RStan), Python (PyStan), or Julia (CmdStan.jl/Stan.jl/StanSample.jl).

## R
To install rstan we first need a C++ compiler. 
To install xcode, open a terminal and issue

`xcode-select --install`

You will get an error if xcode is already installed, in which case you can update xcode via

`softwareupdate --install -a`

Once this is done we install rstan at the R command line via

`install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies=TRUE)`

We will try to fit a simple linear regression to see if Stan is installed correctly. 
Open R and issue:
````
library(rstan)

mod <- "data {
  int N;                  // number of observations
  vector[N] y;            // dependent variable
  vector[N] x;            // independent variable 
}
parameters {
  real beta0;              // intercept
  real beta1;              // slope
  real<lower=1E-15> sigma; // standard deviation of errors
}
model {
  y ~ normal(beta0 + beta1*x, sigma);
}"

N    <- 10
y    <- rnorm(N)
x    <- rnorm(N)
data <- list(y=y,x=x,N=N)

fit  <- stan(model_code=mod, data=data) 

print(fit)
````
This should return you a table of estimated coefficients. 

## Python

We can install pystan directly at terminal via

`pip install pystan`

We will try to fit a simple linear regression to see if Stan is installed correctly. 
Open Python and issue:
````
import numpy
import pystan

mod = """
data {
  int N;                  // number of observations
  vector[N] y;            // dependent variable
  vector[N] x;            // independent variable 
}
parameters {
  real beta0;              // intercept
  real beta1;              // slope
  real<lower=1E-15> sigma; // standard deviation of errors
}
model {
  y ~ normal(beta0 + beta1*x, sigma);
}
"""

N    = 100
y    = numpy.random.normal(0,1,N)
x    = numpy.random.normal(0,1,N)
data = {'y':y, 
        'x':x,
        'N':N}
mod_compile = pystan.StanModel(model_code=mod)
fit = mod_compile.sampling(data=data)
print(fit)
````
This should return you a table of estimated coefficients. 

## Julia

`Stan.jl` uses `cmdstan` which will need to install first.

To do so, open a terminal window, go to the folder where you would like to download `cmdstan`, and clone the `cmdstan` repository as follows.

`git clone https://github.com/stan-dev/cmdstan.git --recursive`

This will download the folder `cmdstan` into whichever folder you currently are.
Now enter that folder and compile `cmdstan`.

```
cd cmdstan
make build
```

This can take 5 or 10 minutes. 

Before starting `julia`, store the path to `cmdstan` in suitable environment variables. 

```
export JULIA_CMDSTAN_HOME=$PWD
export CMDSTAN_HOME=$PWD
```

Once this is done, start `julia` and add `Stan.jl` + `StanSample.jl`.

```
using Pkg
Pkg.add("Stan")
Pkg.add("StanSample")
```

To see if Stan is installed correctly, try to fit a simple linear regression.

````
using StanSample

mod = "
data {
  int N;                  // number of observations
  vector[N] y;            // dependent variable
  vector[N] x;            // independent variable 
}
parameters {
  real beta0;              // intercept
  real beta1;              // slope
  real<lower=1E-15> sigma; // standard deviation of errors
}
model {
  y ~ normal(beta0 + beta1*x, sigma);
}"

data = Dict("N" => 100, "y" => randn(100), "x" => randn(100))
mod_compile = SampleModel("mod", mod, method=StanSample.Sample(save_warmup=true, num_warmup=1000, num_samples=1000, thin=1))
stan_sample(mod_compile, data=data);
println(read_summary(mod_compile));
````

This should return you a table of estimated coefficients. 


# Instructions to install Jupyter notebooks
*NOTE: If you encounter any difficulty with these instructions, please create an issue in GitHub. 
We can help you through the installation if you get stuck and would also like to hear about issues even if you fixed them yourself.*

Issue the following shell command to install jupyter notebooks

```
pip3 install notebook
```
If `pip3` is not installed, see instructions [below](#python-and-pystan). 

# Stan installation instructions

## R and Rstan

### R for Ubuntu (for more general Linux instructions, see below)

If it is not already installed, install R. The easiest way to do so is to use the official Ubuntu repositories (which is often not the latest version of R):
```
sudo apt install r-recommended build-essential
```
For a more up to date version of R follow the instructions [here](https://linuxize.com/post/how-to-install-r-on-ubuntu-18-04/).

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

## Linux
### R for general Linux

Here are R installation instructions for [CentOS](https://linuxize.com/post/how-to-install-r-on-centos-7/) and [Debian](https://linuxize.com/post/how-to-install-r-on-debian-9/), general Linux instructions can be found [here](https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-R-be-installed-_0028Unix_002dlike_0029). 

### RStan

To install Rstan, open the R console as root
```
sudo -i R
```
and then follow the official [installation instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

## Python and pystan

### Python for Ubuntu (for more general Linux instructions, see below)

If it is not already installed, install Python 3 and pip (see below) using the shell command:

```
sudo apt install python3 python3-pip
```
Now, the required python packages can be installed using pip, the Python package installer. In Ubuntu, `pip3` is the Python 3 version of pip. 


### Python for general Linux

[This article](https://www.tecmint.com/install-pip-in-linux/) describes the installation of pip for a variety of Linux distributions. If a Python 3 installation from source is required (typically it is not), instructions can be found [here](https://solarianprogrammer.com/2017/06/30/building-python-ubuntu-wsl-debian/)

### pystan

To install `pystan` (to run Stan), run the follwing shell command to install them:
```
pip3 install pystan
```

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


# Julia
## Julia and JuliaStan

### Julia for Ubuntu (for more general Linux instructions, see below)

While an official Ubuntu repository exists, its version is outdated and we encountered problems running Stan in older versions of Julia. We recommend to follow the general Linux instructions below to get the latest Julia version.

### Julia for general Linux 

To download the latest Linux binary follow the instructions [on the official Julia page](https://julialang.org/downloads/platform.html) (scroll down to the section "Linux and FreeBSD").

### StanSample

#### preparation

To use Stan in Julia, `cmdstan` (a command line version of Stan) is required. It needs to be installed first before downloading any Julia packages. First, find a suitable directory for `cmdstan`, in this example `/some/path/` is used. In the shell execute the following commands (substituting `/some/path/` for the path you selected and created):

```
cd /some/path
git clone https://github.com/stan-dev/cmdstan.git --recursive
cd cmdstan
make build
export JULIA_CMDSTAN_HOME=/some/path/cmdstan
```
The last step above creates a shell variable telling Julia about the location of the `cmdstan` application. Now we are ready to start Julia.

#### in Julia

If not already done, set the `JULIA_CMDSTAN_HOME` variable before starting Julia. In Julia, install the `StanSample` package:

```
import Pkg
Pkg.add("StanSample")
```

--- 
NOTE:

If the above fails with an error message `ERROR: LoadError: No deps.jl file could be found. Please try running Pkg.build("Arpack").` follow the advice and run:

```
Pkg.build("Arpack")
```
Now try again and run:
```
Pkg.add("StanSample")
```
---

We will try to fit a simple linear regression to see if Stan is installed correctly. 
Open Julia and issue:
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

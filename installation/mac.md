# Instructions to install Jupyter notebooks
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

## Python

We can install pystan directly at terminal via

`pip install pystan`


## Julia

Stan.jl uses `cmdstan` which will need to install first.
To do so, clone the following repo from the folder where you would like to install the program

`git clone https://github.com/stan-dev/cmdstan.git --recursive`

which will download the folder `cmdstan` into whichever folder you currently are.
Now navigate into the that folder with

`cd cmdstan`

and issue 

`make build`

which will install `cmdstan`. 
This can take 5 or 10 minutes.

After that installation is finished you can install `Stan.jl`.
To do so, open the Julia prompt and issue

`Pkg.add("Stan")`


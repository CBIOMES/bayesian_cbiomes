# Instructions to install Jupyter notebooks
The first step is to install Anaconda via the instructions here: https://docs.anaconda.com/anaconda/install/.
You can download Anaconda 2 or Anaconda 3.
You will get Jupyter notebooks along with that download.
You will also get a program called 'Anaconda Prompt' which is how we recommend interacting with the notebooks. 
To open or create a notebook, you will open an Anaconda Prompt, navigate to the folder where the notebooks live or are to be created, then type `jupyter notebook` to launch the program.

### Kernels
By default, Anaconda will come with Python kernel pre-loaded in jupyter notebooks so no need to do anything if you are using Python.
You will have a Python 2 kernel if you downloaded Anaconda 2 and a Python 3 kernel if you downloaded Anaconda 3. 

#### R
You can install R at the Anaconda prompt via

`conda install -c r r`

With R installed, you can allow Jupyter notebooks to see R by installing the following  

`conda install -c r r-irkernel`

(You will first need to have R installed on your system). 
Install the R kernel that allows Jupyter notebooks to interact with R from the Anaconda prompt via 

`conda install -c r r-irkernel`

Similarly, you can install any additional package (not contained in the base r-essentials package) via

`conda install -c r <package>`

so to install `rstan` you issue

`conda install -c r rstan`

#### Julia Kernel
Install Julia via: https://julialang.org/downloads/

Start Julia and at the command prompt you can issue the commands 

`using Pkg`
`Pkg.add("IJulia")`

You can close Julia now and launch a Jupyter notebook to see the Julia kernel available. 


# Instructions to install Stan
Instructions to install Stan in Windows will depend on whether you are using R (package: RStan), Python (package: PyStan), or Julia (CmdStan.jl).

### rstan
At the Anaconda prompt issue
`conda install -c r rstan`

You will also have to make sure that you have a C++ compiler that R can interact with. 
You can test this in the R terminal via

`pkgbuild::has_build_tools(debug = TRUE)`

If this returns `FALSE` then you will need to install or update Rtools. 
See instructions here: https://github.com/stan-dev/rstan/wiki/Installing-RStan-from-source-on-Windows

### pystan
You will have to make sure you have a C++ compiler that R can interact with.
At the Anaconda prompt, issue

`conda install libpython m2w64-toolchain -c msys2`

If this fails, see instructions here: https://pystan.readthedocs.io/en/latest/windows.html

Next we will install the dependencies `numpy` and `cython` via

`conda install numpy cython -c conda-forge`

Now we install PyStan via

`pip install python`

or via conda using

`conda install pystan -c conda-forge`

### CmdStan.jl
NOTE: These instructions have failed on some machines for as yet unknown reasons. 
Please contact Greg if you have problems.

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


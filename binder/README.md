Code to allow building working code either in [__Binder__](https://mybinder.org) or on local machine using
[```repo2docker```](https://repo2docker.readthedocs.io/en/latest/). 

To launch repo code running in the free public binder platform use, for example, https://mybinder.org/v2/gh/jpmattern/bayesian_cbiomes/master
This can be somewhat slow, depending on the load from other project on the web that may be using free [project Binder](https://jupyter.org/binder) https://mybinder.org resources at the same time. The virtual machines that are provided for free by https://mybinder.org are not very large. Some of the Python examples in the repo fail due to lack of resources. Interestingly the Julia equivalents work OK.

For heavier cloud use it can be better to set up a dedicated Binder hub resource. The steps for this are described in the [Binderhub documentation](https://binderhub.readthedocs.io/en/latest/zero-to-binderhub/index.html). They are a little involved and require a cloud provider account, but they do work. Most universities have Google, AWS and Azure cloud credit programs. MIT researchers can request Google and Azure credits [here](https://cloud.mit.edu/credits).

Usefully, the binder directory configuration can also be used on a laptop/desktop computer to launch a Docker 
container running the same environment as served by Binder, but executing on a local machine. This can be useful for 
testing things and for executing in an isolated environment. To use on a local computer Docker must be installed on the computer (see https://www.docker.com/get-started) first. Then some somewhat obscure commands can be used to create a 
nicely isolated conda environment for running ```repo2docker``` as follows:

First set up environment in some directory (these commands are for MacOS, but similar commands
will work on Windows or Linux - see https://docs.conda.io/en/latest/miniconda.html )
```
curl  https://repo.continuum.io/miniconda/Miniconda2-4.7.12.1-MacOSX-x86_64.sh > Miniconda2-4.7.12.1-MacOSX-x86_64.sh
chmod +x Miniconda2-4.7.12.1-MacOSX-x86_64.sh
./Miniconda2-4.7.12.1-MacOSX-x86_64.sh -b -p `pwd`/miniconda3
export PATH="`pwd`/miniconda3/bin:$PATH"
. miniconda3/etc/profile.d/conda.sh 
conda create --name myr2d  python=3.6
conda activate myr2d
conda install -c conda-forge jupyter-repo2docker

```

then build and launch via repo2docker 
```
export PATH="`pwd`/miniconda3/bin:$PATH"
. miniconda3/etc/profile.d/conda.sh
conda activate myr2d
repo2docker --user-id 1000 --user-name jovyan https://github.com/jpmattern/bayesian_cbiomes
```

once ```repo2docker``` has finished building the a Docker image, it will launch a container running that image and
text of the form 
```
    To access the notebook, open this file in a browser:
        file:///home/jovyan/.local/share/jupyter/runtime/nbserver-1-open.html
    Or copy and paste one of these URLs:
        http://127.0.0.1:51875/?token=db1d7cfc01cfe66f639b377d918c4121585677ec18415c52
```
will be printed in the terminal running ```repo2docker```. Pasting the URL into a browser on the local
machine should connect to the Jupyter environment for executing the notebooks.

Note - the R dependency in ``runtime.txt`` will be ignored by repo2docker if a Julia config is
requested (specified in REQUIRES). repo2docker can only build one of R or Julia at one time.

__List of files__

- ```REQUIRE``` specifies version of Julia to used.
- ```postBuild``` set of commands that configure Stan for use by Julia and import needed Julia packages. 
- ```requirements.txt``` Python packages that are needed.
- ```runtime.txt``` version of R to install. R will not be installed unless the Julia ```REQUIRE``` file is removed. This is a feature of ```repo2docker```.



__Other Resources__

There are also a variety of other resources available for running code from Notebooks on free cloud resources. These include

- https://www.kaggle.com/
- https://colab.research.google.com/
- https://notebooks.azure.com/
- https://cocalc.com/
- https://datalore.io/
- https://codeocean.com

all of these support Python without any extra steps. Support for R and Julia is available on most of these too, but
generally involves more setup effort and Google searching for tips and tricks!

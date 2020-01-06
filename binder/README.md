Code to allow building working code either in binder or on local machine using
repo2docker. 

To use repo with free public binder use 

To use on a local system 

First set up environment in some directory
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
repo2docker --user-id 1000 --user-name jovyan https://github.com/christophernhill/bayesian_cbiomes

```

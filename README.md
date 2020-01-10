# Bayesian CBIOMES workshop

This repository contains notebooks and files for the Bayesian CBIOMES workshop.

## Preparation

In preparation for the meeting and if Stan is not yet installed on your system follow the [installation instructions](installation/) for your operating system of choice.

## Schedule

![alt text](https://github.com/jpmattern/bayesian_cbiomes/blob/master/schedule.png)

## Tips and tricks

 * Turn any jupyter notebook into a script using
 ```
 jupyter nbconvert --to script notebook_name.ipynb
 ```

 * Some introductions to using Git 
 
   [MITgcm manual Git intro section](https://mitgcm.readthedocs.io/en/latest/contributing/contributing.html#detailed-guide-for-those-less-familiar-with-git-and-github)
   
   [Github guides collection](https://guides.github.com)
   
 * Steps to bring your github master branch up to date with main repo
 
   1. If you have not done so already, download your fork to some directory (e.g. ```~/projects/gits/github.com/YOUR-GITHUB-ORG-NAME/bayesian_cbiomes```) using command such as
   
      ``` 
        $ git clone https://github.com/YOUR-GITHUB-ORG-NAME/bayesian_cbiomes.git 
      ```
      
   1. If you have not done so already, set the _upstream_ repo for your download of the fork
   
      ``` 
        $  cd ~/projects/gits/github.com/YOUR-GITHUB-ORG-NAME/bayesian_cbiomes
        $  git remote add upstream https://github.com/jpmattern/bayesian_cbiomes
      ```
   
   
## For fun
   
 * Some timely videos just added on this years viralest math site on the web [3blue1brown](https://www.3blue1brown.com) recently!
 
   [Intro to Bayes](https://www.youtube.com/watch?v=HZGCoVF3YvM), includes references to the [SS Central America](https://en.wikipedia.org/wiki/SS_Central_America) and other Bayes legends.
   
   [Testing Bayes](https://www.youtube.com/watch?v=U_85TaXbeIo), simple visual working through of Bayes algebra - with nice graphics!



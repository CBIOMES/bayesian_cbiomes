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
   
 * Steps to bring your github fork master branch up to date with main repo
 
   1. If you have not done so already, download your fork to some directory (e.g. ```~/projects/gits/github.com/YOUR-GITHUB-ORG-NAME/bayesian_cbiomes```) using command such as
   
      ``` 
        $ git clone https://github.com/YOUR-GITHUB-ORG-NAME/bayesian_cbiomes.git 
      ```
      
   1. If you have not done so already, set the _upstream_ repo for your download of the fork
   
      ``` 
        $  cd ~/projects/gits/github.com/YOUR-GITHUB-ORG-NAME/bayesian_cbiomes
        $  git remote add upstream https://github.com/jpmattern/bayesian_cbiomes
      ```
      
   1. Now you are set to bring the local copy up to date. 
   
      1. Make sure you are in the relevant directory and have the master branched checked out
   
        ``` 
          $  cd ~/projects/gits/github.com/YOUR-GITHUB-ORG-NAME/bayesian_cbiomes
          $  git checkout master
        ```
        
      2. Now fetch the changes from the upstream repo and merge
        ```
          $ git fetch upstream master
          $ git merge upstream/master
        ```
        
      3. Finally push changes back to your github.com organization
        ```
          $ git push
        ```
        
   1. Awesome, you should now be done! It is possible that the ```git merge upstream/master``` step
      will report a merge conflict error. This happens when the git algorithms can't figure out how
      to merge some changes you have made yourself in the master branch. One practice to manage this
      is to make your own changes on a _feature branch_. The _feature branch_ workflow is described 
      [here](https://www.atlassian.com/git/tutorials/comparing-workflows/feature-branch-workflow).
   
## For fun
   
 * Some timely videos just added on this years viralest math site on the web [3blue1brown](https://www.3blue1brown.com) recently!
 
   [Intro to Bayes](https://www.youtube.com/watch?v=HZGCoVF3YvM), includes references to the [SS Central America](https://en.wikipedia.org/wiki/SS_Central_America) and other Bayes legends.
   
   [Testing Bayes](https://www.youtube.com/watch?v=U_85TaXbeIo), simple visual working through of Bayes algebra - with nice graphics!



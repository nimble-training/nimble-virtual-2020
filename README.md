# nimble-virtual-2020

Materials for the virtual NIMBLE workshop, June 3-5 2020.

To prepare for the workshop:

 - Install NIMBLE (see below)
 - Install additional packages (see below and check back for updates)
 - Join our Slack workspace (highly recommended but not required; see information in email)
 - Download these materials (check back before the workshop on Wednesday for updates)
 - See email for Zoom invitation.

All materials for the workshop are in this Github repository. If you're familiar with Git/Github, you already know how to get all the materials on your computer. If you're not, simply click [here](https://github.com/nimble-training/nimble-virtual-2020/archive/master.zip).

There is some overview information [here](https://htmlpreview.github.io/?https://github.com/nimble-training/nimble-virtual-2020/blob/master/overview.html), including links to the content modules in order.

Time: 9 am - 2 pm California time (GMT-8); noon - 5 pm Eastern US time.
Location: Zoom (see email) and [Slack](https://2020nimbleworkshop.slack.com).

## Slack and Zoom instructions and tips/tricks

Please see [this document](https://docs.google.com/document/d/1hhm6Eco0KevM30aDGdpo0n-gYDcT7IPbw7_vTrb-BVQ/edit?usp=sharing) for discussion of how to use Zoom and Slack during the workshop.

Zoom recordings will be available here (UNDER CONSTRUCTION).

## Schedule (more exact times are in preparation)

Day 1 (Wednesday):

1. Introduction to NIMBLE: Basic concepts and workflows
2. Working with NIMBLE models and converting from WinBUGS/JAGS
3. Comparing and customizing MCMC methods in NIMBLE
4. Strategies for improving MCMC

Day 2 (Thursday):

5. Writing your own functions and distributions 
6. Writing your own MCMC samplers 
7. Spatial modeling (part 1)

Day 3 (Friday):

7. Spatial modeling (part 2)
8. Model selection and Bayesian nonparametrics
9. Special topics breakouts
   - (a) Ecological models
   - (b) Sequential Monte Carlo and maximum likelihood methods

## Help with NIMBLE

Our user manual is [here](https://r-nimble.org/html_manual/cha-welcome-nimble.html).

We have a new 'cheatsheet' [here](https://r-nimble.org/documentation).
(If you have feedback, please let us know during the workshop or by email to nimble.stats@gmail.com)

## Installing NIMBLE

NIMBLE is an R package on CRAN, so in general it will be straightforward to install as with any R package, but you do need a compiler and related tools on your system.  

In summary, here are the steps.

1. Install compiler tools on your system. [https://r-nimble.org/download](https://r-nimble.org/download) has more details on how to install *Rtools* on Windows and how to install the command line tools of *Xcode* on a Mac. Note that if you have packages requiring a compiler (e.g., *Rcpp*) on your computer, you should already have the compiler tools installed.

2. Install the *nimble* package from CRAN in the usual fashion for an R package. More details (including troubleshooting tips) can also be found in Section 4 of the [NIMBLE manual](https://r-nimble.org/html_manual/cha-installing-nimble.html).

3) To test that things are working please run the following code  in R:

```
library(nimble)
code <- nimbleCode({
  y ~ dnorm(0,1)
})
model <- nimbleModel(code)
cModel <- compileNimble(model)
```


If that runs without error, you're all set. If not, please see the troubleshooting tips and email nimble.stats@gmail.com directly if you can't get things going.  

In general we encourage you to update to the most recent version of NIMBLE, 0.9.1.

IMPORTANT: for those of you using R 4.0 under Windows, you MUST use version 0.9.1.

#### (Not required) Development version(s) of NIMBLE

Sometimes we make an update or new feature available on a github branch before it is released.  In the event a need arises to install from a branch, you can do so as follows (for branch "devel"):

```
library(remotes)
install_github("nimble-dev/nimble", ref = "devel", subdir = "packages/nimble")
```

## Installing additional packages

Some of the packages we will use (beyond those automatically installed with `nimble`) can be installed as follows:

```
install.packages(c("mcmcplots", "CARBayesdata", "sp", "spdep", "classInt"))
```

`compareMCMCs` is a package in development that is not yet on CRAN:

```
library(remotes)
install_github("nimble-dev/compareMCMCs", subdir = "compareMCMCs")
```

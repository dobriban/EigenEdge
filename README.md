# EigenEdge
Computing with Eigenvalue Distributions of Large Random Matrices (including *Spectrode*)

The **EigenEdge** MATLAB package contains open source implementations 
of methods for working with eigenvalue distributions of large random matrices 
of the covariance type. 


## Contents 

* a detailed documentation with examples (`\Readme`)
* methods to compute the limit empirical spectrum of sample covariance matrices: **Spectrode**, fixed point and Newton. *Spectrode* is described in: 
*Dobriban, E.* : Efficient Computation of Limit Spectra of Sample Covariance Matrices, submitted.  http://arxiv.org/abs/1504.02935
* methods to compute arbitrary moments and quantiles of the limit spectrum 
* Matlab scripts to reproduce all computational results of the above paper (see `\Experiments\Spectrode` folder)

An altenative to installing the package is to download only stand-alone scripts: 
* For *Spectrode*, the appropriate file is at `\Code\compute_esd_ode`. Please see the readme for a usage example. 

This package is work in progress. Suggestions and comments are welcome.

# EigenEdge
Computing with Eigenvalue Distributions of Large Random Matrices

The **EigenEdge** MATLAB package contains open source implementations 
of methods for working with eigenvalue distributions of large random matrices 
of the covariance type (known as general Marchenko-Pastur distributions). 


## Contents 

* a documentation with examples: `\Readme\readme.pdf`
* the SPECTRODE method to compute the limit empirical spectrum of sample covariance matrices (equivalently: compute the Marchenko-Pastur forward map) 
* methods to compute moments and quantiles of the limit spectrum 
* optimal linear spectral statistics (LSS) for testing in Principal Component Analysis

## Additional contents

This package also has the software to reproduce the computational results of the following papers of the author: 
* Dobriban. E,  *Efficient Computation of Limit Spectra of Sample Covariance Matrices*, Random Matrices: Theory Appl., 04, 1550019 (2015). http://arxiv.org/abs/1507.01649
* Dobriban. E,  *Sharp detection in PCA under correlations: all eigenvalues matter*, http://arxiv.org/abs/1602.06896

## Installation

For full notes on installing the package, please see the readme. 

##Additional notes
This package is work in progress. Suggestions and comments are welcome.
An R version is under development, and is available from https://github.com/dobriban/Spectrode-R/

#Example
A typical example eigenvalue distribution computed with this software looks as follows: 
![Alt text](https://github.com/dobriban/EigenEdge/blob/master/Experiments/Examples/Illustration_mixture_3.png?raw=true "Optional Title")
This example is explained in detail in the readme at `\Readme\readme.pdf`. The input is the distribution of population eigenvalues, which is a mixture of point masses and a uniform density. The output is the distribution of sample eigenvalues, which has a smooth density on several disjoint intervals. There is a close match with a Monte Carlo simulation.


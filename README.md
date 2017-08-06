# EigenEdge
The **EigenEdge** MATLAB package contains open source implementations 
of methods for working with eigenvalue distributions of large random matrices. 
The focus is on covariance-type, or "general Marchenko-Pastur" distributions. 

These tools can be used in high-dimensional statistics, wireless communications, and finance, among other areas.

If you are interested in trying this out, have suggestions, or think this may have applications to your field, please do not hesitate to get in touch with me. 

## Contents 

* a documentation with examples: [readme.pdf](https://github.com/dobriban/EigenEdge/blob/master/readme/readme.pdf)
* the SPECTRODE method to compute the limit empirical spectrum of sample covariance matrices (equivalently: compute the Marchenko-Pastur forward map) 
* methods to compute moments and quantiles of the limit spectrum 
* optimal linear spectral statistics (LSS) for testing in Principal Component Analysis
* spiked models: "forward" and "inverse" spike maps in standard and general spiked models

## References

This package is partially based on methods proposed in the following papers: 
* Dobriban. E,  *Efficient Computation of Limit Spectra of Sample Covariance Matrices*, Random Matrices: Theory Appl., 04, 1550019 (2015). [pdf](https://github.com/dobriban/Papers/blob/master/Dobriban%20-%20Efficient%20computation%20of%20limit%20spectra%20of%20sample%20covariance%20matrices-2015-RMTA.pdf)
* Dobriban. E,  *Sharp detection in PCA under correlations: all eigenvalues matter*, The Annals of Statistics. [pdf](https://github.com/dobriban/Papers/blob/master/Dobriban%20-%20Sharp%20detection%20in%20PCA%20under%20correlations.%20all%20eigenvalues%20matter%20-%202017%20-%20AoS.pdf)

The package also has the software to reproduce the computational results of the above papers. 

## Installation

For full notes on installing the package, please see the readme. 

##Additional notes
This package is work in progress. Suggestions and comments are welcome.
An R version is under development, and is available from https://github.com/dobriban/Spectrode-R/

##Acknowledgements
The author is very grateful to the following indivuals for their invaluable advice, help, comments, bug submissions, etc. with this package:
* Romain Couillet, Iain Johnstone, Jack Silverstein.
* The students in STATS 325 "Random Matrix Theory and High-Dimensional Statistics" at Stanford University, Spring 2016. In particular: Jeha Yang, 	Simon Rosenberg. 

#Example
A typical example eigenvalue distribution computed with this software looks as follows: 
![Alt text](https://github.com/dobriban/EigenEdge/blob/master/Experiments/Examples/Illustration_mixture_3.png?raw=true "Optional Title")
This example is explained in detail in the readme at `\Readme\readme.pdf`. The input is the distribution of population eigenvalues, which is a mixture of point masses and a uniform density. The output is the distribution of sample eigenvalues, which has a smooth density on several disjoint intervals. There is a close match with a Monte Carlo simulation.


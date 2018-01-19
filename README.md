[![Travis-CI Build Status](https://travis-ci.org/SMAC-Group/ns-gmwm.svg?branch=master)](https://travis-ci.org/SMAC-Group/ns-gmwm)

# ns-gmwm

This package is in development and is made available for experimental purposes.
It is a minor extension of the `gmwm` package. It is used for the paper 
``Estimation of Inertial Sensor Stochastic Characteristics under Varying 
Environmental Conditions'' by S. Orso, P. Clausen, S. Guerrier, J. Skaloud.

## Install Instructions

First, install the R specific package dependencies, then the package itself.
``` bash
# Install dependencies
install.packages(c("RcppArmadillo","devtools"))

# Install the package
devtools::install_github("SMAC-Group/ns-gmwm")
```
Note that other requirements might be necessary depending on your system (more information in the future).

We recommend to install the `gmwm` package.
``` bash
devtools::install_github("SMAC-Group/gmwm")
```

## Usage
In its actual version, the package has only one functionnality: compute the following objective function   
<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{1}{2}\frac{1}{K}\sum_{k=1}^K\lVert\hat{\nu}_k&space;-&space;\nu(\theta,x_k)\rVert^2_{\Omega}" target="_blank"><img src="https://latex.codecogs.com/png.latex?\frac{1}{2}\frac{1}{K}\sum_{k=1}^K\lVert\hat{\nu}_k&space;-&space;\nu(\theta,x_k)\rVert^2_{\Omega}" title="\frac{1}{2}\frac{1}{K}\sum_{k=1}^K\lVert\hat{\nu}_k - \nu(\theta,x_k)\rVert^2_{\Omega}" /></a>   

and its gradient     

<a href="https://www.codecogs.com/eqnedit.php?latex=\frac{1}{K}\sum_{k=1}^K\frac{\partial}{\partial\theta^T}\nu(\theta,x_k)\Omega[\hat{\nu}_k-\nu(\theta,x_k)]" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\frac{1}{K}\sum_{k=1}^K\frac{\partial}{\partial\theta^T}\nu(\theta,x_k)\Omega[\hat{\nu}_k-\nu(\theta,x_k)]" title="\frac{1}{K}\sum_{k=1}^K\frac{\partial}{\partial\theta^T}\nu(\theta,x_k)\Omega[\hat{\nu}_k-\nu(\theta,x_k)]" /></a>

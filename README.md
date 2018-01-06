# ns-gmwm

This package is in development and is made available for experimental purposes.
It is a minor extension of the `gmwm` package.

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
In its actual version, the package has only one functionnality: compute the following objective function $$\frac{1}{2}\frac{1}{K}\sum_{k=1}^K\lVert\hat{\nu}_k - \nu(\theta,x_k)\rVert^2_{\Omega}} $$ and its gradient $$ \frac{1}{K}\sum_{k=1}^K(\hat{\nu}_k - \nu(\theta,x_k))\Omega\frac{\partial}{\partial\theta}\nu(\theta,x_k)$$


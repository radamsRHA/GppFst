# GppFst: Genomic Posterior Predictive Fst simulations  
**NOTE See the file <https://github.com/radamsRHA/GppFst/blob/master/vignettes/GppFst_Tutorial.pdf> for detailed instructions**

Steps for conducting posterior predictive simulations for Fst and Dxy. 

1. Estimate joint posterior distributions for coalescent parameters of a two population coalescent model (population size parameters 4*Ne*u, divergence times). This can be accomplished using a number of programs that estimate mutispecies coalescent model parameters.
2. Input a sample of this posterior distribution using the function `read.posterior`.
3. Input the empirical locus coverage (for both populations) and locus length for all loci
4. Conduct Fst and Dxy PPS using the input MCMC and empirical pararmeters  

## Installing R package GppFst from github

The R package GppFst is freely available to download and distribute from github <https://github.com/radamsRHA/GppFst/>. To install and load GppFst, you must first install the R packages `devtools`, `phybase`, and `Geneland`.

```
install.packages("devtools")
install.packages("Geneland")
```
The packages, `phybase` should be downloaded <http://faculty.franklin.uga.edu/lliu/content/phybase>.  Save it to your desktop and then use this code to install it from source:

```
install.packages("~/Desktop/phybase", repos = NULL, type = "source")
```

Now using devtools we can install `GppFst` from github:

```
library(devtools)
install_github("radamsRHA/GppFst")
library(GppFst) # Load package GppFst
library(phybase) # Load dependancy phybase
library(Geneland) # Load dependancy Geneland
```

UPDATE: it seems that Geneland has been removed from Cran. As a potential work around, you can install gfortran for mac, for example (https://github.com/fxcoudert/gfortran-for-macOS/releases), and then use the command:

```
install_github("https://github.com/cran/Geneland")
```

To begin using `GppFst` try using our vignette with example files provided with this package. See the directory ./GppFst_Tutorial/ for example files and tutorial. 



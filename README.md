[![DOI](https://zenodo.org/badge/119578865.svg)](https://zenodo.org/doi/10.5281/zenodo.12706763)

# scrattch
### **S**ingle-**c**ell **R**NA-seq **a**nalysis for **t**ranscriptomic **t**ype **ch**aracterization

This is the umbrella package for the `scrattch` suite of R packages from the Allen Institute for Brain Science. It is modeled after the [`tidyverse`](https://www.tidyverse.org/) package.  

You can use `scrattch` to automatically install or update these packages:  

[`scrattch.io`](https://github.com/AllenInstitute/scrattch.io) - File handling and data formatting  
[`scrattch.hicat`](https://github.com/AllenInstitute/scrattch.hicat) - **H**ierarchical, **i**terative **c**lustering for **a**nalysis of  **t**ranscriptomics  
[`scrattch.bigcat`](https://github.com/AllenInstitute/scrattch.bigcat) - Clustering analysis for extremely large single cell dataset  
[`scrattch.mapping`](https://github.com/AllenInstitute/scrattch.mapping) - Generalized mapping scripts for RNA-seq, Patch-seq or any gene expression data  
[`scrattch.vis`](https://github.com/AllenInstitute/scrattch.vis) - Plotting functions for visualization of RNA-seq data  
[`tasic2016data`](https://github.com/AllenInstitute/tasic2016data) - Data from Tasic, et al. (2016), which is used for demos  

If you're interested in only one of these modules, you can install them separately. `scrattch` will install them all together.  

### Installation

To install all scrattch packages along with their Github and BioConductor dependencies, use:
```
devtools::install_github("AllenInstitute/scrattch")

scrattch::install_scrattch()
```

The previous, internal version of `scrattch` has been split to two packages:  
File handling and data formatting are now part of [`scrattch.io`](https://github.com/AllenInstitute/scrattch.io).  
Plotting and visualization are now part of [`scrattch.vis`](https://github.com/AllenInstitute/scrattch.vis).  

Should you need the previous version, it can still be installed using:  
```
devtools::install_github("AllenInstitute/scrattch", ref = "archive")
```

### Function List

There are now a lot of functions available in `scrattch` packages. To assist in finding which scrattch package a function is stored in, you can check this CSV file stored in the `scrattch` umbrella package:  
[scrattch function list](https://github.com/AllenInstitute/scrattch/blob/dev/inst/scrattch_function_list.csv)

### License

The license for this package is available on Github at: [https://github.com/AllenInstitute/scrattch/blob/master/LICENSE](https://github.com/AllenInstitute/scrattch/blob/master/LICENSE)

### Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/scrattch/blob/master/CONTRIBUTION

### Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.


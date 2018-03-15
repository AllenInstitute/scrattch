# scrattch
### **S**ingle-**c**ell **R**NA-seq **a**nalysis for **t**ranscriptomic **t**ype **ch**aracterization

### **Note: This is currently a placeholder to outline expected functionality, and not yet fully functional!**

This is the umbrella package for the `scrattch` suite of R packages from the Allen Institute for Brain Science. It is modeled after the [`tidyverse`](https://www.tidyverse.org/) package.  

You can use `scrattch` to automatically install or update these packages:  

[`scrattch.io`](https://github.com/AllenInstitute/scrattch.io) - File handling and data formatting  
[`scrattch.iterclust`](https://github.com/AllenInstitute/iterclust) - Iterative clustering for cell type characterization  
[`scrattch.lowcat`](https://github.com/AllenInstitute/lowcat)  - Low-coverage accessibility and transcriptomics data analysis (scATAC-seq and scRNA-seq comparisons)  
[`scrattch.vis`](https://github.com/AllenInstitute/scrattch.vis) - Plotting and data visualization  
[`tasic2016data`](https://github.com/AllenInstitute/tasic2016data) - Data from Tasic, et al. (2016), which is used for demos  

If you're interested in only one of these modules, you can install them separately. `scrattch` will install them all together.  

To install all scrattch packages and their CRAN dependencies, use:
```
devtools::install_github("AllenInstitute/scrattch")
install_scrattch()
```

Once installed, you can retrieve all Bioconductor dependencies using:
```
scrattch::install_bioc_deps()
```

The previous internal version of `scrattch` has been split to two packages:  
File handling and data formatting are now part of [`scrattch.io`](https://github.com/AllenInstitute/scrattch.io).  
Plotting and visualization are now part of [`scrattch.vis`](https://github.com/AllenInstitute/scrattch.vis).  


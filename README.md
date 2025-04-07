[![DOI](https://zenodo.org/badge/119578865.svg)](https://zenodo.org/doi/10.5281/zenodo.12706763)

# scrattch
## **S**ingle-**c**ell **R**NA-seq **a**nalysis for **t**ranscriptomic **t**ype **ch**aracterization

This is the umbrella package for the `scrattch` suite of R packages from the Allen Institute for Brain Science. It is modeled after the [`tidyverse`](https://www.tidyverse.org/) package.  You can use `scrattch` to automatically install or update some of the underlying packages and can run the remaining packages in docker environments.

<img src="https://github.com/user-attachments/assets/6c29a501-6934-486f-8b8e-6b72a21a9b6c" width="200" />

## Scrattch packages

Scrattch includes several packages for clustering, mapping, and data formatting and visualization, along with example data for demos.  These include:

**Data preparation:** file formats and schema

* [`scrattch.taxonomy`](https://github.com/AllenInstitute/scrattch.taxonomy) - Taxonomy building scripts for RNA-seq based taxonomies following the [Allen Institute (AIT) schema](https://github.com/AllenInstitute/AllenInstituteTaxonomy/tree/main/schema).  
* [`scrattch.io`](https://github.com/AllenInstitute/scrattch.io) - [deprecated]. Library for file handling and data formatting, replaced by `scrattch.taxonomy` in 2024  

**Data analysis:** cell clustering and mapping (also called label transfer)

* [`scrattch.hicat`](https://github.com/AllenInstitute/scrattch.hicat) - **H**ierarchical, **i**terative **c**lustering for **a**nalysis of  **t**ranscriptomics  
* [`scrattch.bigcat`](https://github.com/AllenInstitute/scrattch.bigcat) - Clustering analysis for extremely large single cell dataset  
* [`scrattch.mapping`](https://github.com/AllenInstitute/scrattch.mapping) - Generalized mapping scripts for single cell RNA-seq, Patch-seq, spatial transcriptomics, or related data types  
* [`scrattch.patchseq`](https://github.com/AllenInstitute/scrattch.patchseq) - Functions for generating additional QC metrics and output files for patch-seq analysis   

**Data visualization**

* [`scrattch.vis`](https://github.com/AllenInstitute/scrattch.vis) - Plotting functions for visualization of single cell RNA-seq data  

**Example data:** small RNA-seq data sets

* [`tasic2016data`](https://github.com/AllenInstitute/tasic2016data) - Data from [Tasic, et al. (2016)](https://pubmed.ncbi.nlm.nih.gov/26727548/), which is used for demos  
* [`hodge2019data`](https://github.com/AllenInstitute/hodge2019data) - Data subset from [Hodge, et al. (2019)](https://pubmed.ncbi.nlm.nih.gov/31435019/), which is used for demos  

If you're interested in only one of these modules, you can install them separately. That said, we recommend using the installation instructions below to install combinations of `scrattch` packages to ensure they interact properly.  

### Related content

Several related websites and R and python libraries are outside of the `scrattch` suite, but are either used as part of `scrattch` libraries or directly work with `scrattch` outputs.  These include (but are not limited to):

* [`bmark`](https://github.com/AllenInstitute/bmark) - Standardized strategies for benchmarking clustering and mapping results  
* [`transcriptomic_clustering`](https://github.com/AllenInstitute/transcriptomic_clustering) - Python implementation of scrattch.hicat clustering  
* [`cell_type_mapper`](https://github.com/AllenInstitute/cell_type_mapper) - Python implementation of hierarchical mapping algorithm used in `scrattch.mapping` and [`MapMyCells`](https://portal.brain-map.org/atlases-and-data/bkp/mapmycells)  
* [`ACE`](https://github.com/AllenInstitute/ace) - R Shiny and [web-based](https://sea-ad.shinyapps.io/ACEapp/) app for comparison of annotations, including clustering and mapping results  
* [`mfishtools`](https://github.com/AllenInstitute/mfishtools) - Functions for gene selection and analysis of spatial transcriptomics data  


## Installation

**We strongly encourage the use of docker to the `scrattch` suite.** In particular, several functions in `scrattch.taxonomy` and `scrattch.mapping` have known issues in certain R environments.  That said, we provide options for installing and running R in both a docker environment and through standard R approaches.

### Using docker (RECOMMENDED)

The current docker version is accessible through [Docker Hub](https://hub.docker.com/u/alleninstitute). As of 26 March 2025 the [Docker version](https://hub.docker.com/r/alleninst/scrattch/tags) is **`docker://alleninst/scrattch:1.1.2`**.  This corresponds to AIT (v1.1.2) (see [the Allen Institute Taxonomy GitHub respository](https://github.com/AllenInstitute/AllenInstituteTaxonomy/) for details).

Docker can be run on some HPC environments that use singularity as follows:

* **Non-interactive**: `singularity shell --cleanenv [Docker version] Rscript YOUR_CODE.R`
* **Interactive**: `singularity shell --cleanenv [Docker version]`
* **To create a sif file for use in other environments**: `singularity pull scrattch:[#.#.#].sif [Docker version]`

If you cannot figure out how to use Docker in your specific environment, please post an issue.

*--WARNING-- The 1.1.2 docker listed above provides all the tooling for AIT (v1.1.2) and some functionality for scrattch.mapping, but is broken for hierarchical mapping and all scrattch.patchseq functionality. An update mid-April will bring both of these packages back up to speed with the AIT schema / format.* 

### Running `scrattch` in R

While we advise using the provided docker, you can install all `scrattch` packages along with their GitHub and BioConductor dependencies, as follows:
```
devtools::install_github("AllenInstitute/scrattch")
scrattch::install_scrattch()
```
Note that `doMC` may need to be installed manually from the download link at https://r-forge.r-project.org/R/?group_id=947 if you use Windows.

### Installing previous versions

Two historical versions of scrattch are included in this package. These can be safely run without using docker, but are missing several recent components of the `scrattch` suite.

* **scrattch_2023** is the stable version of the package prior to the release of `scrattch.mapping`, `scrattch.taxonomy`, `scrattch.patchseq`, and `hodge2019data`  
* **archive** is the original package from ~2018, and should not be used for most folks

Should you need one of these previous versions, they can still be installed using:  
```
devtools::install_github("AllenInstitute/scrattch", ref = "scrattch_2023") # -OR-
devtools::install_github("AllenInstitute/scrattch", ref = "archive")
```


## Usage and contribution

`Scrattch` is under active development. Please reach out if you have any challenges or suggestions!

### Documentation

There is no specific documentation for `scrattch`. However, each child package has extensive documentation available via help commands within R and/or through the corresponding GitHub page.  Conversion of documentation for select packages to ReadTheDocs in planned for April.

### License

The license for this package is available on Github at: [https://github.com/AllenInstitute/scrattch/blob/master/LICENSE](https://github.com/AllenInstitute/scrattch/blob/master/LICENSE)

### Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/scrattch/blob/master/CONTRIBUTION

### Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.  We encourage community involvement in child packages directly, rather than through the `scrattch` umbrella package, when appropriate.


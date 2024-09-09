
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ScReNI: single-cell regulatory network inference through integrating scRNA-seq and scATAC-seq data

ScReNI is a novel algorithm that leverages k-nearest neighbors and
random forest algorithms to integrate single-cell RNA sequencing
(scRNA-seq) and single-cell ATAC sequencing (scATAC-seq) data, enabling
the inference of gene regulatory networks at the single-cell level.
Designed to be versatile, ScReNI is adept at handling both paired and
unpaired datasets from scRNA-seq and scATAC-seq. It stands out with its
capability of identifying cell-enriched regulators based on each
cell-specific network.

ScReNI has the following five key steps:

-   **Cells measured by scRNA-seq and scATAC-seq are clustered
    together**

-   **Neighbors for each cell are identified using the k-nearest
    neighbors algorithm**

-   **The weights of transcription factors regulating target genes are
    calculated using random forest**

-   **ScReNI is validated by assessing cell-specific regulatory
    networks**

-   **Cell-enriched regulators are identified based on each
    cell-specific network**

<img src="Readme%20figure/ScReNI_schematics.png"
style="width:70.0%;height:70.0%" />

## News

### Version1.0.0: First version of ScReNI was published!

### Installation

Install ScReNI from Github:

``` r
install.packages('devtools')
devtools::install_github('Xuxl2020/ScReNI')
```

or download source file and install locally

``` r
devtools::install_github('path/to/ScReNI.tar.gz')
```

To check if ScReNI is installed correctly, run:

``` r
library(ScReNI)
```

## Example

We use unpaired scRNA-seq and scATAC-seq data from retinal development
as example datasets.

[Using ScReNI to analyze unpaired scRNA-seq and scATAC-seq data from
retinal
development](https://htmlpreview.github.io/?https://github.com/Xuxl2020/ScReNI/blob/master/docs/ScReNI_tutorial.html)

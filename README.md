# GeoPAT 2.0

## Overview

GeoPAT 2.0 (Geospatial Pattern Analysis Toolbox) is a standalone suite of modules written in C and dedicated to analysis of large Earth Science datasets in their entirety using spatial and/or temporal patterns. 
Global scale, high resolution spatial datasets are available but are mostly used in small pieces for local studies. 
GeoPAT enables studying them in their entirety.
GeoPAT’s core idea is to tessellate global spatial data into grid of square blocks of original cells (pixels).
This transforms data from its original form (huge number of cells each having simple content) to a new form (much smaller number of supercells/blocks with complex content).
Complex cell contains a pattern of original variable.
GeoPAT provides means for succinct description of such patterns and for calculation of similarity between patterns.
This enables spatial analysis such as search, change detection, segmentation, and clustering to be performed on the grid of complex cells (local patterns).

![](https://github.com/Nowosad/geopat2_manual/raw/master/figs/logo.png)

The GeoPAT 2.0 software is available at http://sil.uc.edu/cms/index.php?id=geopat2.

The GeoPAT 2.0 manual is at https://rawgit.com/Nowosad/geopat2_manual/master/output/GeoPAT2_Manual.pdf.

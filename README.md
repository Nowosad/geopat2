# GeoPAT 2

[![Build Status](https://travis-ci.org/Nowosad/geopat2.svg?branch=master)](https://travis-ci.org/Nowosad/geopat2)

## Overview

**GeoPAT** 2 (Geospatial Pattern Analysis Toolbox) is a standalone suite of modules written in C and dedicated to analysis of large Earth Science datasets in their entirety using spatial and/or temporal patterns. 
Global scale, high resolution spatial datasets are available but are mostly used in small pieces for local studies. 
**GeoPAT** enables studying them in their entirety.
**GeoPAT**’s core idea is to tessellate global spatial data into grid of square blocks of original cells (pixels).
This transforms data from its original form (huge number of cells each having simple content) to a new form (much smaller number of supercells/blocks with complex content).
Complex cell contains a pattern of original variable.
**GeoPAT** provides means for succinct description of such patterns and for calculation of similarity between patterns.
This enables spatial analysis such as search, change detection, segmentation, and clustering to be performed on the grid of complex cells (local patterns).

![](https://github.com/Nowosad/geopat2_manual/raw/master/figs/logo.png)

The **GeoPAT** 2 software and manual are available at http://sil.home.amu.edu.pl/index.php?id=geopat2.

## Installation

Installation instruction in detail can be found in the [GeoPAT 2 manual](https://rawgit.com/Nowosad/geopat2_manual/master/output/GeoPAT2_Manual.pdf). 

### Windows

The installer for Windows x64 is available at https://github.com/Nowosad/geopat2win/raw/master/GPAT2setup.exe

### MacOS

The GeoPAT 2 Unix executable programs, compiled for Mac OS with the Apple M1 Max processor (an ARM64 architecture), are available at https://github.com/Nowosad/geopat2mac
These programs can be run from a shell script, or with R's `system` command.

### Building from source code

To build **GeoPAT** 2 from the source code, the development files for GDAL are required.
They can be installed on Ubuntu using:

```bash
sudo add-apt-repository ppa:ubuntugis/ubuntugis-unstable --yes
sudo apt-get --yes --force-yes update -qq
sudo apt-get install --yes libgdal-dev
```
... or on Fedora with:

```bash
sudo dnf install gdal-devel
```

The source code of **GeoPAT** 2 is available at https://github.com/Nowosad/geopat2/archive/master.zip. 
This archive should be unpacked, compiled and installed, e.g. with:

```bash
wget https://github.com/Nowosad/geopat2/archive/master.zip
unzip master.zip
mv geopat2-master geopat2
cd geopat2
make
sudo make install
```

## How to cite **GeoPAT** 2

**Netzel P., Nowosad J., Jasiewicz J., Niesterowicz J, Stepinski T.**, 2018. *GeoPAT 2: user's manual*. Zenodo. http://doi.org/10.5281/zenodo.1291123

## Blog posts

1. [GeoPAT 2: Software for Pattern-Based Spatial and Temporal Analysis](http://nowosad.github.io/post/geopat-2-software-for-pattern-based-spatial-and-temporal-analysis)
2. [Pattern-based Spatial Analysis - core ideas](http://nowosad.github.io/post/pattern-based-spatial-analysis-core-ideas)
3. [Finding similar local landscapes](http://nowosad.github.io/post/geopat-2-search)
4. [Quantifying temporal change of landscape pattern](http://nowosad.github.io/post/geopat-2-compare)
5. [Pattern-based regionalization](http://nowosad.github.io/post/geopat-2-segmentation)
6. [Moving beyond pattern-based analysis: Additional applications of GeoPAT 2](http://nowosad.github.io/post/geopat-2-extend)
7. [GeoPAT2: Entropy calculations for local landscapes](https://nowosad.github.io/post/geopat-2-ent/)

## Workshop at GEOSTAT 2018

- [Video](https://www.youtube.com/watch?v=_yKbVqR_Zfc)
- [Slides](https://nowosad.github.io/geostat18/geostat18_nowosad)

<!--
## **GeoPAT** related papers

List of the papers related to the **GeoPAT** software. The preprints can be found at http://sil.home.amu.edu.pl/index.php?id=journal-papers.

- Nowosad, J. and Stepinski, T. F., 2018. Towards machine ecoregionalization of Earth's landmass using pattern segmentation method. International Journal of Applied Earth Observation and Geoinformation, 69, pp.110-118.
- Nowosad J., Stepinski  T. F., 2018. Global inventory of landscape patterns and latent variables of landscape spatial configuration, Ecological Indicators 89, pp. 159-167.
- Niesterowicz, J. and Stepinski, T.F., 2017. Pattern-based, multi-scale segmentation and regionalization of EOSD land cover. International Journal of Applied Earth Observation and Geoinformation, 62, pp.192-200.
- Niesterowicz, J., Stepinski, T.F. and Jasiewicz, J., 2016. Unsupervised regionalization of the United States into landscape pattern types. International Journal of Geographical Information Science, 30(7), pp.1450-1468.
- Netzel, P. and Stepinski, T.F., 2017. World Climate Search and Classification Using a Dynamic Time Warping Similarity Function. In Advances in Geocomputation (pp. 181-195). Springer
- Netzel, P. and Stepinski, T., 2016. On using a clustering approach for global climate classification. Journal of Climate, 29(9), pp.3387-3401.
- Niesterowicz, J., Stepinski, T. and Jasiewicz, J., 2016, January. Unsupervised Delineation of Urban Structure Types Using High Resolution RGB Imagery. In International Conference on GIScience Short Paper Proceedings (Vol. 1, No. 1).
- Jasiewicz, J., Niesterowicz, J. and Stepinski, T., 2016, January. Multi-resolution, pattern-based segmentation of very large raster datasets. In International Conference on GIScience Short Paper Proceedings (Vol. 1, No. 1).
- Netzel, P., Jasiewicz, J. and Stepinski, T., 2016, January. TerraEx–a GeoWeb app for world-wide content-based search and distribution of elevation and landforms data. In International Conference on GIScience Short Paper Proceedings (Vol. 1, No. 1).
- Netzel, P. and Stepinski, T.F., 2015. Pattern-Based Assessment of Land Cover Change on Continental Scale With Application to NLCD 2001–2006. IEEE Transactions on Geoscience and Remote Sensing, 53(4), pp.1773-1781.
- Jasiewicz, J., Netzel, P. and Stepinski, T., 2015. GeoPAT: A toolbox for pattern-based information retrieval from large geospatial databases. Computers & Geosciences, 80, pp.62-73.
- Stepinski, T.F., Niesterowicz, J. and Jasiewicz, J., 2015. Pattern-based regionalization of large geospatial datasets using complex object-based image analysis. Procedia Computer Science, 51, pp.2168-2177.
- Jasiewicz, J., Netzel, P. and Stepinski, T.F., 2014. Landscape similarity, retrieval, and machine mapping of physiographic units. Geomorphology, 221, pp.104-112.
- Jasiewicz, J., Netzel, P. and Stepinski, T.F., 2014, July. Retrieval of pattern-based information from giga-cells categorical rasters—Concept and new software. In Geoscience and Remote Sensing Symposium (IGARSS), 2014 IEEE International (pp. 1785-1788). IEEE.
- Netzel, P. and Stepinski, T.F., 2014, July. Pattern-based assessment of 2001/2006 land cover change over the entire United States. In Geoscience and Remote Sensing Symposium (IGARSS), 2014 IEEE International (pp. 4188-4191). IEEE.
- Stepinski, T.F., Netzel, P. and Jasiewicz, J., 2014. LandEx—a GeoWeb tool for query and retrieval of spatial patterns in land cover datasets. IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 7(1), pp.257-266.
- Jasiewicz, J. and Stepinski, T.F., 2013. Example-based retrieval of alike land-cover scenes from NLCD2006 database. IEEE Geoscience and Remote Sensing Letters, 10(1), pp.155-159.
-->

## Acknowledgments

This work was supported by the University of Cincinnati Space Exploration Institute and by the grant NNX15AJ47G from the National Aeronautics and Space Administration (NASA).
We also want to thank D G Rossiter for preparing and sharing the MacOS version of **GeoPAT** 2.

Data Processing
================

``` r
library(raster)
library(tidyverse)
```

General details
---------------

We are creating a dataset of social-ecological covariates at a range of spatial resolutions: 1km, 5km, 10km, 25km, 50km, 100km. All data are open.

Land-cover data
---------------

We use [Corine Land Cover 2012](https://land.copernicus.eu/pan-european/corine-land-cover/clc-2012) at 100m resolution to calculate land-cover diversity at each resolution.

Elevation data
--------------

We use [EU-DEM v1.1](https://land.copernicus.eu/pan-european/satellite-derived-products/eu-dem/eu-dem-v1.1/view) at 30m resolution to calculate topographic variability at each resolution.

Population data
---------------

We use [Population density disaggregated with Corine land cover 2000](https://www.eea.europa.eu/data-and-maps/data/population-density-disaggregated-with-corine-land-cover-2000-2) at xx resolution to calculate total population at each resolution.
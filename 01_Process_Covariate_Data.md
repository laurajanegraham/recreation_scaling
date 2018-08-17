Process Covariate Data
================

``` r
# Libraries
library(raster)
library(tidyverse)
library(gdalUtils)
library(doParallel)
library(knitr)

# Analysis resolutions
rln <- c(1, 5, 10, 25, 50, 100)

# Calculate the number of cores
no_cores <- detectCores() - 1
```

General details
---------------

We are creating a dataset of social-ecological covariates at a range of spatial resolutions: 1km, 5km, 10km, 25km, 50km, 100km. All data are open.

Elevation data
--------------

We use [EU-DEM v1.1](https://land.copernicus.eu/pan-european/satellite-derived-products/eu-dem/eu-dem-v1.1/view) at 100m resolution to calculate topographic variability at each resolution. Note that we are using a version that we previously aggregated (mean) to 100m resolution to reduce processing time. We may wish to revisit.

``` r
dem <- raster("~/DATA/PHYSICAL/elev/eu_dem_1.1/eudem_100m_aggregated.tif")

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, c("dem"))
clusterEvalQ(cl, "raster")

# aggregate to each resolution
parLapply(cl, rln, function(x) {
  fname <- paste0("data/dem_", x, "km.tif")
  raster::aggregate(dem, (x*1000)/100, 
            fun=var, 
            filename = fname,
            overwrite = TRUE)
})

stopCluster(cl)
```

Land-cover data
---------------

We use [Corine Land Cover 2012](https://land.copernicus.eu/pan-european/corine-land-cover/clc-2012) at 100m resolution to calculate land-cover diversity at each resolution.

``` r
align_rasters(unaligned = "C:/Users/lg1u16/DATA/LULC/clc2012/g100_clc12_V18_5.tif", 
              reference = "C:/Users/lg1u16/DATA/PHYSICAL/elev/eu_dem_1.1/eudem_100m_aggregated.tif",
              dstfile = "data/clc_100m.tif",
              r = "mode",
              nThreads = "ALL_CPUS",
              overwrite = TRUE)

clc <- raster("data/clc_100m.tif")

 
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, c("clc"))
clusterEvalQ(cl, c("raster", "grainchanger"))

# aggregate to each resolution
parLapply(cl, rln, function(x) {
  fname <- paste0("data/clc_", x, "km.tif")
  raster::aggregate(clc, (x*1000)/100, 
            fun=function(y, ...) grainchanger::diversity(y, lc_class = 1:44), 
            filename = fname,
            overwrite = TRUE)
})

stopCluster(cl)
```

Population data
---------------

We use [Population density disaggregated with Corine land cover 2000](https://www.eea.europa.eu/data-and-maps/data/population-density-disaggregated-with-corine-land-cover-2000-2) at 100m resolution to calculate total population at each resolution.

``` r
align_rasters(unaligned = "C:/Users/lg1u16/DATA/SOCIAL/EEA_gridded_population/popu01clcv5.tif", 
                     reference = "C:/Users/lg1u16/DATA/PHYSICAL/elev/eu_dem_1.1/eudem_100m_aggregated.tif",
                     dstfile = "data/pop_100m.tif",
                     nThreads = "ALL_CPUS")
```

    ## NULL

``` r
pop <- raster("data/pop_100m.tif")

# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, c("pop"))
clusterEvalQ(cl, "raster")
```

    ## [[1]]
    ## [1] "raster"
    ## 
    ## [[2]]
    ## [1] "raster"
    ## 
    ## [[3]]
    ## [1] "raster"

``` r
# aggregate to each resolution
parLapply(cl, rln, function(x) {
  fname <- paste0("data/pop_", x, "km.tif")
  raster::aggregate(pop, (x*1000)/100, 
            fun=sum, 
            filename = fname,
            overwrite = TRUE)
})
```

    ## [[1]]
    ## class       : RasterLayer 
    ## dimensions  : 3878, 3878, 15038884  (nrow, ncol, ncell)
    ## resolution  : 1000, 1000  (x, y)
    ## extent      : 2603525, 6481525, 1494900, 5372900  (xmin, xmax, ymin, ymax)
    ## coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
    ## data source : C:\Users\lg1u16\GIT_PROJECTS\SCALEFORES\recreation_scaling\data\pop_1km.tif 
    ## names       : pop_1km 
    ## values      : 0, 4072700  (min, max)
    ## 
    ## 
    ## [[2]]
    ## class       : RasterLayer 
    ## dimensions  : 776, 776, 602176  (nrow, ncol, ncell)
    ## resolution  : 5000, 5000  (x, y)
    ## extent      : 2603525, 6483525, 1492900, 5372900  (xmin, xmax, ymin, ymax)
    ## coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
    ## data source : C:\Users\lg1u16\GIT_PROJECTS\SCALEFORES\recreation_scaling\data\pop_5km.tif 
    ## names       : pop_5km 
    ## values      : 0, 66105183  (min, max)
    ## 
    ## 
    ## [[3]]
    ## class       : RasterLayer 
    ## dimensions  : 388, 388, 150544  (nrow, ncol, ncell)
    ## resolution  : 10000, 10000  (x, y)
    ## extent      : 2603525, 6483525, 1492900, 5372900  (xmin, xmax, ymin, ymax)
    ## coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
    ## data source : C:\Users\lg1u16\GIT_PROJECTS\SCALEFORES\recreation_scaling\data\pop_10km.tif 
    ## names       : pop_10km 
    ## values      : 0, 187955505  (min, max)
    ## 
    ## 
    ## [[4]]
    ## class       : RasterLayer 
    ## dimensions  : 156, 156, 24336  (nrow, ncol, ncell)
    ## resolution  : 25000, 25000  (x, y)
    ## extent      : 2603525, 6503525, 1472900, 5372900  (xmin, xmax, ymin, ymax)
    ## coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
    ## data source : C:\Users\lg1u16\GIT_PROJECTS\SCALEFORES\recreation_scaling\data\pop_25km.tif 
    ## names       : pop_25km 
    ## values      : 0, 530242916  (min, max)
    ## 
    ## 
    ## [[5]]
    ## class       : RasterLayer 
    ## dimensions  : 78, 78, 6084  (nrow, ncol, ncell)
    ## resolution  : 50000, 50000  (x, y)
    ## extent      : 2603525, 6503525, 1472900, 5372900  (xmin, xmax, ymin, ymax)
    ## coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
    ## data source : C:\Users\lg1u16\GIT_PROJECTS\SCALEFORES\recreation_scaling\data\pop_50km.tif 
    ## names       : pop_50km 
    ## values      : 0, 755303241  (min, max)
    ## 
    ## 
    ## [[6]]
    ## class       : RasterLayer 
    ## dimensions  : 39, 39, 1521  (nrow, ncol, ncell)
    ## resolution  : 1e+05, 1e+05  (x, y)
    ## extent      : 2603525, 6503525, 1472900, 5372900  (xmin, xmax, ymin, ymax)
    ## coord. ref. : +proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs 
    ## data source : C:\Users\lg1u16\GIT_PROJECTS\SCALEFORES\recreation_scaling\data\pop_100km.tif 
    ## names       : pop_100km 
    ## values      : 0, 1046067791  (min, max)

``` r
stopCluster(cl)
```

Session Info
------------

``` r
session <- devtools::session_info()
session[[1]]
```

    ##  setting  value                       
    ##  version  R version 3.5.0 (2018-04-23)
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United Kingdom.1252 
    ##  tz       Europe/London               
    ##  date     2018-08-16

``` r
session[[2]] %>% kable
```

| package     | \*  | version  | date       | source                             |
|:------------|:----|:---------|:-----------|:-----------------------------------|
| assertthat  |     | 0.2.0    | 2017-04-11 | CRAN (R 3.5.0)                     |
| backports   |     | 1.1.2    | 2017-12-13 | CRAN (R 3.5.0)                     |
| base        | \*  | 3.5.0    | 2018-04-23 | local                              |
| bindr       |     | 0.1.1    | 2018-03-13 | CRAN (R 3.5.0)                     |
| bindrcpp    |     | 0.2.2    | 2018-03-29 | CRAN (R 3.5.0)                     |
| broom       |     | 0.4.4    | 2018-03-29 | CRAN (R 3.5.0)                     |
| cellranger  |     | 1.1.0    | 2016-07-27 | CRAN (R 3.5.0)                     |
| cli         |     | 1.0.0    | 2017-11-05 | CRAN (R 3.5.0)                     |
| codetools   |     | 0.2-15   | 2016-10-05 | CRAN (R 3.5.0)                     |
| colorspace  |     | 1.3-2    | 2016-12-14 | CRAN (R 3.5.0)                     |
| compiler    |     | 3.5.0    | 2018-04-23 | local                              |
| crayon      |     | 1.3.4    | 2017-09-16 | CRAN (R 3.5.0)                     |
| datasets    | \*  | 3.5.0    | 2018-04-23 | local                              |
| devtools    |     | 1.13.5   | 2018-02-18 | CRAN (R 3.5.0)                     |
| digest      |     | 0.6.15   | 2018-01-28 | CRAN (R 3.5.0)                     |
| doParallel  | \*  | 1.0.11   | 2017-09-28 | CRAN (R 3.5.0)                     |
| dplyr       | \*  | 0.7.6    | 2018-06-29 | CRAN (R 3.5.1)                     |
| evaluate    |     | 0.10.1   | 2017-06-24 | CRAN (R 3.5.0)                     |
| forcats     | \*  | 0.3.0    | 2018-02-19 | CRAN (R 3.5.0)                     |
| foreach     | \*  | 1.4.4    | 2017-12-12 | CRAN (R 3.5.0)                     |
| foreign     |     | 0.8-70   | 2017-11-28 | CRAN (R 3.5.0)                     |
| gdalUtils   | \*  | 2.0.1.14 | 2018-04-23 | CRAN (R 3.5.0)                     |
| ggplot2     | \*  | 3.0.0    | 2018-07-03 | CRAN (R 3.5.1)                     |
| glue        |     | 1.3.0    | 2018-07-17 | CRAN (R 3.5.1)                     |
| graphics    | \*  | 3.5.0    | 2018-04-23 | local                              |
| grDevices   | \*  | 3.5.0    | 2018-04-23 | local                              |
| grid        |     | 3.5.0    | 2018-04-23 | local                              |
| gtable      |     | 0.2.0    | 2016-02-26 | CRAN (R 3.5.0)                     |
| haven       |     | 1.1.1    | 2018-01-18 | CRAN (R 3.5.0)                     |
| hms         |     | 0.4.2    | 2018-03-10 | CRAN (R 3.5.0)                     |
| htmltools   |     | 0.3.6    | 2017-04-28 | CRAN (R 3.5.0)                     |
| httr        |     | 1.3.1    | 2017-08-20 | CRAN (R 3.5.0)                     |
| iterators   | \*  | 1.0.9    | 2017-12-12 | CRAN (R 3.5.0)                     |
| jsonlite    |     | 1.5      | 2017-06-01 | CRAN (R 3.5.0)                     |
| knitr       | \*  | 1.20     | 2018-02-20 | CRAN (R 3.5.0)                     |
| lattice     |     | 0.20-35  | 2017-03-25 | CRAN (R 3.5.0)                     |
| lazyeval    |     | 0.2.1    | 2017-10-29 | CRAN (R 3.5.0)                     |
| lubridate   |     | 1.7.4    | 2018-04-11 | CRAN (R 3.5.0)                     |
| magrittr    |     | 1.5      | 2014-11-22 | CRAN (R 3.5.0)                     |
| memoise     |     | 1.1.0    | 2017-04-21 | CRAN (R 3.5.0)                     |
| methods     | \*  | 3.5.0    | 2018-04-23 | local                              |
| mnormt      |     | 1.5-5    | 2016-10-15 | CRAN (R 3.5.0)                     |
| modelr      |     | 0.1.2    | 2018-05-11 | CRAN (R 3.5.0)                     |
| munsell     |     | 0.4.3    | 2016-02-13 | CRAN (R 3.5.0)                     |
| nlme        |     | 3.1-137  | 2018-04-07 | CRAN (R 3.5.0)                     |
| parallel    | \*  | 3.5.0    | 2018-04-23 | local                              |
| pillar      |     | 1.2.2    | 2018-04-26 | CRAN (R 3.5.0)                     |
| pkgconfig   |     | 2.0.1    | 2017-03-21 | CRAN (R 3.5.0)                     |
| plyr        |     | 1.8.4    | 2016-06-08 | CRAN (R 3.5.0)                     |
| psych       |     | 1.8.4    | 2018-05-06 | CRAN (R 3.5.0)                     |
| purrr       | \*  | 0.2.5    | 2018-05-29 | CRAN (R 3.5.1)                     |
| R.methodsS3 |     | 1.7.1    | 2016-02-16 | CRAN (R 3.5.0)                     |
| R.oo        |     | 1.22.0   | 2018-04-22 | CRAN (R 3.5.0)                     |
| R.utils     |     | 2.6.0    | 2017-11-05 | CRAN (R 3.5.0)                     |
| R6          |     | 2.2.2    | 2017-06-17 | CRAN (R 3.5.0)                     |
| raster      | \*  | 2.6-7    | 2017-11-13 | CRAN (R 3.5.1)                     |
| Rcpp        |     | 0.12.18  | 2018-07-23 | CRAN (R 3.5.1)                     |
| readr       | \*  | 1.1.1    | 2017-05-16 | CRAN (R 3.5.0)                     |
| readxl      |     | 1.1.0    | 2018-04-20 | CRAN (R 3.5.0)                     |
| reshape2    |     | 1.4.3    | 2017-12-11 | CRAN (R 3.5.0)                     |
| rgdal       |     | 1.3-4    | 2018-08-03 | CRAN (R 3.5.1)                     |
| rlang       |     | 0.2.1    | 2018-05-30 | CRAN (R 3.5.1)                     |
| rmarkdown   |     | 1.10     | 2018-06-11 | CRAN (R 3.5.0)                     |
| rprojroot   |     | 1.3-2    | 2018-01-03 | CRAN (R 3.5.0)                     |
| rstudioapi  |     | 0.7      | 2017-09-07 | CRAN (R 3.5.0)                     |
| rvest       |     | 0.3.2    | 2016-06-17 | CRAN (R 3.5.0)                     |
| scales      |     | 0.5.0    | 2017-08-24 | CRAN (R 3.5.0)                     |
| sp          | \*  | 1.2-7    | 2018-01-19 | CRAN (R 3.5.0)                     |
| stats       | \*  | 3.5.0    | 2018-04-23 | local                              |
| stringi     |     | 1.1.7    | 2018-03-12 | CRAN (R 3.5.0)                     |
| stringr     | \*  | 1.3.1    | 2018-05-10 | CRAN (R 3.5.0)                     |
| tibble      | \*  | 1.4.2    | 2018-01-22 | CRAN (R 3.5.0)                     |
| tidyr       | \*  | 0.8.0    | 2018-01-29 | CRAN (R 3.5.0)                     |
| tidyselect  |     | 0.2.4    | 2018-02-26 | CRAN (R 3.5.0)                     |
| tidyverse   | \*  | 1.2.1    | 2017-11-14 | CRAN (R 3.5.0)                     |
| tools       |     | 3.5.0    | 2018-04-23 | local                              |
| utils       | \*  | 3.5.0    | 2018-04-23 | local                              |
| withr       |     | 2.1.2    | 2018-07-20 | Github (<jimhester/withr@fe56f20>) |
| xml2        |     | 1.2.0    | 2018-01-24 | CRAN (R 3.5.0)                     |
| yaml        |     | 2.1.19   | 2018-05-01 | CRAN (R 3.5.0)                     |

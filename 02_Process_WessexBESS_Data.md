Process Wessex BESS data
================

``` r
# Libraries
library(raster)
library(tidyverse)
library(gdalUtils)
library(doParallel)
library(knitr)
library(sf)

# create study site for the WBESS data
# currently based on a 10km resolution. Change if we change the resolutions at all
study_ext <- c(xmin = 3463525, xmax = 3523525, ymin = 3162900, ymax = 3232900)
```

Wessex BESS recreation visits data
----------------------------------

These data were collected using a PPGIS methodology as part of the Wessex BESS project. NB there are more details available (e.g. type of visit). The bounding box for the Wessex BESS area is roughly 75 x 60 km, so the maximum resolution we will use is 25km (although even this is too coarse).

``` r
visits_sf <- st_read("data/response/All_PPGIS_points_LR.shp", quiet = TRUE) %>% 
  st_transform(crs = 3035)

visits_sp <- as(visits_sf, "Spatial")
```

We need to convert from points to raster at the 6 analysis resolutions (NB we won't be using all of them, this is for consistency across studies for now).

``` r
# loop through the dem files to get rasters for the visits
fnames <- list.files("data/covariates", pattern = "dem_mean", full.names = TRUE)

for(f in fnames) {
  ras <- raster(f)
  ras <- crop(ras, study_ext)
  rln <- res(ras)/1000
  fname <- paste0("data/response/wbess_", rln[1], "km.tif")
  rec_ras <- rasterize(visits_sp, 
                       ras, 
                       'LUCY_FID', 
                       fun=function(x,...)length(x), 
                       filename = fname,
                       overwrite = TRUE)
}
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
    ##  date     2018-09-17

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
| class       |     | 7.3-14   | 2015-08-30 | CRAN (R 3.5.0)                     |
| classInt    |     | 0.2-3    | 2018-04-16 | CRAN (R 3.5.0)                     |
| cli         |     | 1.0.0    | 2017-11-05 | CRAN (R 3.5.0)                     |
| codetools   |     | 0.2-15   | 2016-10-05 | CRAN (R 3.5.0)                     |
| colorspace  |     | 1.3-2    | 2016-12-14 | CRAN (R 3.5.0)                     |
| compiler    |     | 3.5.0    | 2018-04-23 | local                              |
| crayon      |     | 1.3.4    | 2017-09-16 | CRAN (R 3.5.0)                     |
| datasets    | \*  | 3.5.0    | 2018-04-23 | local                              |
| DBI         |     | 1.0.0    | 2018-05-02 | CRAN (R 3.5.0)                     |
| devtools    |     | 1.13.5   | 2018-02-18 | CRAN (R 3.5.0)                     |
| digest      |     | 0.6.15   | 2018-01-28 | CRAN (R 3.5.0)                     |
| doParallel  | \*  | 1.0.11   | 2017-09-28 | CRAN (R 3.5.0)                     |
| dplyr       | \*  | 0.7.6    | 2018-06-29 | CRAN (R 3.5.1)                     |
| e1071       |     | 1.6-8    | 2017-02-02 | CRAN (R 3.5.0)                     |
| evaluate    |     | 0.10.1   | 2017-06-24 | CRAN (R 3.5.0)                     |
| forcats     | \*  | 0.3.0    | 2018-02-19 | CRAN (R 3.5.0)                     |
| foreach     | \*  | 1.4.4    | 2017-12-12 | CRAN (R 3.5.0)                     |
| foreign     |     | 0.8-71   | 2018-07-20 | CRAN (R 3.5.1)                     |
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
| sf          | \*  | 0.6-3    | 2018-05-17 | CRAN (R 3.5.1)                     |
| sp          | \*  | 1.2-7    | 2018-01-19 | CRAN (R 3.5.0)                     |
| spData      |     | 0.2.8.3  | 2018-03-25 | CRAN (R 3.5.0)                     |
| stats       | \*  | 3.5.0    | 2018-04-23 | local                              |
| stringi     |     | 1.1.7    | 2018-03-12 | CRAN (R 3.5.0)                     |
| stringr     | \*  | 1.3.1    | 2018-05-10 | CRAN (R 3.5.0)                     |
| tibble      | \*  | 1.4.2    | 2018-01-22 | CRAN (R 3.5.0)                     |
| tidyr       | \*  | 0.8.0    | 2018-01-29 | CRAN (R 3.5.0)                     |
| tidyselect  |     | 0.2.4    | 2018-02-26 | CRAN (R 3.5.0)                     |
| tidyverse   | \*  | 1.2.1    | 2017-11-14 | CRAN (R 3.5.0)                     |
| tools       |     | 3.5.0    | 2018-04-23 | local                              |
| udunits2    |     | 0.13     | 2016-11-17 | CRAN (R 3.5.0)                     |
| units       |     | 0.5-1    | 2018-01-08 | CRAN (R 3.5.0)                     |
| utils       | \*  | 3.5.0    | 2018-04-23 | local                              |
| withr       |     | 2.1.2    | 2018-07-20 | Github (<jimhester/withr@fe56f20>) |
| xml2        |     | 1.2.0    | 2018-01-24 | CRAN (R 3.5.0)                     |
| yaml        |     | 2.1.19   | 2018-05-01 | CRAN (R 3.5.0)                     |

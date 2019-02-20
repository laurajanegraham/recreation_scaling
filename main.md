Scale-dependency in drivers of outdoor recreation
================

# Data

Response data are MENE visits, and flickr visits.

For each response type and resolution, the covariates are:

  - Proportion of agricultural land cover (clc\_agri; supply)
  - Proporiton of natural land covers (clc\_prop; supply)
  - Diversity of land-cover types (clc\_shei; supply)
  - Amount of protected area (pa; supply)
  - Mean elevation (dem\_mean; supply)
  - Total population (pop; demand)
  - Distance to nearest major town or city (dist; demand)

## Data exploration

### Distribution, sample size and coverage of recreation proxies

![](main_files/figure-gfm/rec_distribution-1.png)<!-- -->

Number of non-zero cells (total is total number of cells with either a
flickr or MENE observation):

| resolution | total | flickr | mene |
| :--------- | ----: | -----: | ---: |
| 5km        |  2193 |   2193 | 2193 |
| 10km       |  1042 |   1042 | 1042 |
| 25km       |   207 |    207 |  207 |
| 50km       |    51 |     51 |   51 |

### Spatial coverage of the observations

![](main_files/figure-gfm/rec_spatial-1.png)<!-- -->

### Distribution of covariates

![](main_files/figure-gfm/cov_dist-1.png)<!-- -->

  - clc\_agri: left skew at all resolutions
  - clc\_shei: fine
  - all others: right skew at all resolutions

We log transform clc\_prop, dem\_range, pa, dist and pop; we square
transform clc\_agri.

### Correlation between variables

    ## [[1]]

![](main_files/figure-gfm/cov_corrs-1.png)<!-- -->

    ## 
    ## [[2]]

![](main_files/figure-gfm/cov_corrs-2.png)<!-- -->

    ## 
    ## [[3]]

![](main_files/figure-gfm/cov_corrs-3.png)<!-- -->

    ## 
    ## [[4]]

![](main_files/figure-gfm/cov_corrs-4.png)<!-- -->

High correlation between variables along with small sample size might
mean dropping the 50km resolution analysis (for now care to be taken
with this analysis).

# Q1: How similar are proxies for recreation at different scales?

Hypothesis: The MENE visits and Flickr photos are proxies for different
kinds of recreation. MENE visits are day-to-day local recreation; Flickr
photos represent \[tourism, aesthetic, holiday etc. - need a name for
this\].

Test: MENE visits will not explain a large proportion of the variation
in Flickr photos,

![](main_files/figure-gfm/comp_proxy-1.png)<!-- -->

We have one rather large outlier at each resolution - remove this and
test:

![](main_files/figure-gfm/comp_proxy_rmoutliers-1.png)<!-- -->

We are fitting a poisson glm for each resolution with flickr photos as
response, and mene visits as covariate. We plot \(D^2\), coefficient
estimate and spatial plot of residuals. We do this for the dataset with
and without the outliers.

![](main_files/figure-gfm/comp_proxy_mod-1.png)<!-- -->

The outliers make a huge difference to the variance explained. We are
currently removing the largest outlier at each resolution in subsequent
analyses. NB the outlier at each resolutions is in London.

Spatial distribution of the residuals

![](main_files/figure-gfm/comp_proxy_spatial-1.png)<!-- -->

When the residuals are positive, the number of flickr photos is much
higher than expected based on the MENE visits. The residuals are more
often positive than negative, and by a bigger margin. This means that
when they diverge from each other, it’s because there are more flickr
photos than expected. We could do a more detailed analysis, but looks
like the higher positive residuals are mainly around National
Parks.

# Q2: What are the relative importances of drivers of supply of and demand for ecosystem services at different scales?

Hypothesis: supply will be a more important predictor of the flickr
photos, demand of the MENE visits.

Test: variance partitioning

## Model selection

We fit 3 negative binomial models for each resolution and response type.
The models are supply, demand and full (combined supply and demand
models):

  - Proportion of agricultural land cover (clc\_agri; supply)
  - Proporiton of natural land covers (clc\_prop; supply)
  - Diversity of land-cover types (clc\_shei; supply)
  - Mean elevation (dem\_\_mean\_range; supply)
  - Proportion of area covered by protected area (pa; supply)
  - Total population (pop; demand)
  - Distance to nearest urban area (dist;
supply)

| resolution | dataset | mod\_full\_aic | mod\_demand\_aic | mod\_supply\_aic |
| :--------- | :------ | -------------: | ---------------: | ---------------: |
| 5km        | mene    |      8957.1226 |        9146.4748 |        9235.0444 |
| 10km       | mene    |      5318.0472 |        5456.8823 |        5704.7983 |
| 25km       | mene    |      1560.2132 |        1607.0113 |        1722.6162 |
| 50km       | mene    |       485.1326 |         492.0830 |         541.9825 |
| 5km        | flickr  |     19708.7308 |       20612.6152 |       19861.8065 |
| 10km       | flickr  |     11092.4701 |       11640.0854 |       11317.6449 |
| 25km       | flickr  |      2798.4752 |        2945.8813 |        2871.6089 |
| 50km       | flickr  |       789.7393 |         833.5414 |         814.9753 |

In all cases the full (supply & demand) model is the best performing
model as judged by AIC. At all resolutions, there is more support for
the demand than supply model when MENE observations are used as a proxy
for recreation. Similarly, at all resolutions there is more support for
the supply than demand model when flickr photos are used as a proxy.

### Variance partitioning

Total variance explained and how it’s partitioned for the full
model

![](main_files/figure-gfm/varpart-1.png)<!-- -->![](main_files/figure-gfm/varpart-2.png)<!-- -->

  - When using Flickr as a proxy, the supply part of the model counts
    for \> 50% of the explained variance at all resolutions.
  - When using MENE visits as a proxy, the supply part of the model
    never accounts for more than -% of the explained variance
  - The full model explains slightly less variance in the flickr data
    than it does the MENE data
  - Scale dependencies: in both cases, the coarser the scale, the more
    variance is explained by demand for rather than supply of the
    service (although this is very marginal in the flickr model).

### Variable relationships

![](main_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

## Session Info

    ##  setting  value                       
    ##  version  R version 3.5.1 (2018-07-02)
    ##  os       Windows 10 x64              
    ##  system   x86_64, mingw32             
    ##  ui       RTerm                       
    ##  language (EN)                        
    ##  collate  English_United Kingdom.1252 
    ##  ctype    English_United Kingdom.1252 
    ##  tz       Europe/London               
    ##  date     2018-12-05

|              | package      | ondiskversion | loadedversion | path                                           | loadedpath                                     | attached | is\_base | date       | source                                 | md5ok | library                            |
| ------------ | :----------- | :------------ | :------------ | :--------------------------------------------- | :--------------------------------------------- | :------- | :------- | :--------- | :------------------------------------- | :---- | :--------------------------------- |
| assertthat   | assertthat   | 0.2.0         | 0.2.0         | C:/Users/lg1u16/R/win-library/3.5/assertthat   | C:/Users/lg1u16/R/win-library/3.5/assertthat   | FALSE    | FALSE    | 2017-04-11 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| backports    | backports    | 1.1.2         | 1.1.2         | C:/Users/lg1u16/R/win-library/3.5/backports    | C:/Users/lg1u16/R/win-library/3.5/backports    | FALSE    | FALSE    | 2017-12-13 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| base64enc    | base64enc    | 0.1.3         | 0.1-3         | C:/Users/lg1u16/R/win-library/3.5/base64enc    | C:/Users/lg1u16/R/win-library/3.5/base64enc    | FALSE    | FALSE    | 2015-07-28 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| bindr        | bindr        | 0.1.1         | 0.1.1         | C:/Users/lg1u16/R/win-library/3.5/bindr        | C:/Users/lg1u16/R/win-library/3.5/bindr        | FALSE    | FALSE    | 2018-03-13 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| bindrcpp     | bindrcpp     | 0.2.2         | 0.2.2         | C:/Users/lg1u16/R/win-library/3.5/bindrcpp     | C:/Users/lg1u16/R/win-library/3.5/bindrcpp     | TRUE     | FALSE    | 2018-03-29 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| broom        | broom        | 0.5.0         | 0.5.0         | C:/Users/lg1u16/R/win-library/3.5/broom        | C:/Users/lg1u16/R/win-library/3.5/broom        | TRUE     | FALSE    | 2018-07-17 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| callr        | callr        | 3.0.0         | 3.0.0         | C:/Users/lg1u16/R/win-library/3.5/callr        | C:/Users/lg1u16/R/win-library/3.5/callr        | FALSE    | FALSE    | 2018-08-24 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| cellranger   | cellranger   | 1.1.0         | 1.1.0         | C:/Users/lg1u16/R/win-library/3.5/cellranger   | C:/Users/lg1u16/R/win-library/3.5/cellranger   | FALSE    | FALSE    | 2016-07-27 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| class        | class        | 7.3.14        | 7.3-14        | C:/Program Files/R/R-3.5.1/library/class       | C:/Program Files/R/R-3.5.1/library/class       | FALSE    | FALSE    | 2015-08-30 | CRAN (R 3.5.1)                         | NA    | C:/Program Files/R/R-3.5.1/library |
| classInt     | classInt     | 0.2.3         | 0.2-3         | C:/Users/lg1u16/R/win-library/3.5/classInt     | C:/Users/lg1u16/R/win-library/3.5/classInt     | FALSE    | FALSE    | 2018-04-16 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| cli          | cli          | 1.0.1         | 1.0.1         | C:/Users/lg1u16/R/win-library/3.5/cli          | C:/Users/lg1u16/R/win-library/3.5/cli          | FALSE    | FALSE    | 2018-09-25 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| codetools    | codetools    | 0.2.15        | 0.2-15        | C:/Program Files/R/R-3.5.1/library/codetools   | C:/Program Files/R/R-3.5.1/library/codetools   | FALSE    | FALSE    | 2016-10-05 | CRAN (R 3.5.1)                         | TRUE  | C:/Program Files/R/R-3.5.1/library |
| colorspace   | colorspace   | 1.3.2         | 1.3-2         | C:/Users/lg1u16/R/win-library/3.5/colorspace   | C:/Users/lg1u16/R/win-library/3.5/colorspace   | FALSE    | FALSE    | 2016-12-14 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| corrr        | corrr        | 0.3.0         | 0.3.0         | C:/Users/lg1u16/R/win-library/3.5/corrr        | C:/Users/lg1u16/R/win-library/3.5/corrr        | TRUE     | FALSE    | 2018-07-31 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| crayon       | crayon       | 1.3.4         | 1.3.4         | C:/Users/lg1u16/R/win-library/3.5/crayon       | C:/Users/lg1u16/R/win-library/3.5/crayon       | FALSE    | FALSE    | 2017-09-16 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| DBI          | DBI          | 1.0.0         | 1.0.0         | C:/Users/lg1u16/R/win-library/3.5/DBI          | C:/Users/lg1u16/R/win-library/3.5/DBI          | FALSE    | FALSE    | 2018-05-02 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| debugme      | debugme      | 1.1.0         | 1.1.0         | C:/Users/lg1u16/R/win-library/3.5/debugme      | C:/Users/lg1u16/R/win-library/3.5/debugme      | FALSE    | FALSE    | 2017-10-22 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| desc         | desc         | 1.2.0         | 1.2.0         | C:/Users/lg1u16/R/win-library/3.5/desc         | C:/Users/lg1u16/R/win-library/3.5/desc         | FALSE    | FALSE    | 2018-05-01 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| devtools     | devtools     | 2.0.0         | 2.0.0         | C:/Users/lg1u16/R/win-library/3.5/devtools     | C:/Users/lg1u16/R/win-library/3.5/devtools     | FALSE    | FALSE    | 2018-10-19 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| digest       | digest       | 0.6.18        | 0.6.18        | C:/Users/lg1u16/R/win-library/3.5/digest       | C:/Users/lg1u16/R/win-library/3.5/digest       | FALSE    | FALSE    | 2018-10-10 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| dplyr        | dplyr        | 0.7.7         | 0.7.7         | C:/Users/lg1u16/R/win-library/3.5/dplyr        | C:/Users/lg1u16/R/win-library/3.5/dplyr        | TRUE     | FALSE    | 2018-10-16 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| e1071        | e1071        | 1.7.0         | 1.7-0         | C:/Users/lg1u16/R/win-library/3.5/e1071        | C:/Users/lg1u16/R/win-library/3.5/e1071        | FALSE    | FALSE    | 2018-07-28 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| evaluate     | evaluate     | 0.12          | 0.12          | C:/Users/lg1u16/R/win-library/3.5/evaluate     | C:/Users/lg1u16/R/win-library/3.5/evaluate     | FALSE    | FALSE    | 2018-10-09 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| forcats      | forcats      | 0.3.0         | 0.3.0         | C:/Users/lg1u16/R/win-library/3.5/forcats      | C:/Users/lg1u16/R/win-library/3.5/forcats      | TRUE     | FALSE    | 2018-02-19 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| fs           | fs           | 1.2.6         | 1.2.6         | C:/Users/lg1u16/R/win-library/3.5/fs           | C:/Users/lg1u16/R/win-library/3.5/fs           | FALSE    | FALSE    | 2018-08-23 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| gamlss       | gamlss       | 5.1.2         | 5.1-2         | C:/Users/lg1u16/R/win-library/3.5/gamlss       | C:/Users/lg1u16/R/win-library/3.5/gamlss       | TRUE     | FALSE    | 2018-10-06 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| gamlss.data  | gamlss.data  | 5.1.0         | 5.1-0         | C:/Users/lg1u16/R/win-library/3.5/gamlss.data  | C:/Users/lg1u16/R/win-library/3.5/gamlss.data  | TRUE     | FALSE    | 2018-10-06 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| gamlss.dist  | gamlss.dist  | 5.1.0         | 5.1-0         | C:/Users/lg1u16/R/win-library/3.5/gamlss.dist  | C:/Users/lg1u16/R/win-library/3.5/gamlss.dist  | TRUE     | FALSE    | 2018-10-06 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| GGally       | GGally       | 1.4.0         | 1.4.0         | C:/Users/lg1u16/R/win-library/3.5/GGally       | C:/Users/lg1u16/R/win-library/3.5/GGally       | TRUE     | FALSE    | 2018-05-17 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| ggplot2      | ggplot2      | 3.1.0         | 3.1.0         | C:/Users/lg1u16/R/win-library/3.5/ggplot2      | C:/Users/lg1u16/R/win-library/3.5/ggplot2      | TRUE     | FALSE    | 2018-10-25 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| glue         | glue         | 1.3.0         | 1.3.0         | C:/Users/lg1u16/R/win-library/3.5/glue         | C:/Users/lg1u16/R/win-library/3.5/glue         | FALSE    | FALSE    | 2018-07-17 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| gtable       | gtable       | 0.2.0         | 0.2.0         | C:/Users/lg1u16/R/win-library/3.5/gtable       | C:/Users/lg1u16/R/win-library/3.5/gtable       | FALSE    | FALSE    | 2016-02-26 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| gtools       | gtools       | 3.8.1         | 3.8.1         | C:/Users/lg1u16/R/win-library/3.5/gtools       | C:/Users/lg1u16/R/win-library/3.5/gtools       | TRUE     | FALSE    | 2018-06-26 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| haven        | haven        | 1.1.2         | 1.1.2         | C:/Users/lg1u16/R/win-library/3.5/haven        | C:/Users/lg1u16/R/win-library/3.5/haven        | FALSE    | FALSE    | 2018-06-27 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| hier.part    | hier.part    | 1.0.4         | 1.0-4         | C:/Users/lg1u16/R/win-library/3.5/hier.part    | C:/Users/lg1u16/R/win-library/3.5/hier.part    | TRUE     | FALSE    | 2013-01-08 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| highr        | highr        | 0.7           | 0.7           | C:/Users/lg1u16/R/win-library/3.5/highr        | C:/Users/lg1u16/R/win-library/3.5/highr        | FALSE    | FALSE    | 2018-06-09 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| hms          | hms          | 0.4.2         | 0.4.2         | C:/Users/lg1u16/R/win-library/3.5/hms          | C:/Users/lg1u16/R/win-library/3.5/hms          | FALSE    | FALSE    | 2018-03-10 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| htmltools    | htmltools    | 0.3.6         | 0.3.6         | C:/Users/lg1u16/R/win-library/3.5/htmltools    | C:/Users/lg1u16/R/win-library/3.5/htmltools    | FALSE    | FALSE    | 2017-04-28 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| httr         | httr         | 1.3.1         | 1.3.1         | C:/Users/lg1u16/R/win-library/3.5/httr         | C:/Users/lg1u16/R/win-library/3.5/httr         | FALSE    | FALSE    | 2017-08-20 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| jsonlite     | jsonlite     | 1.5           | 1.5           | C:/Users/lg1u16/R/win-library/3.5/jsonlite     | C:/Users/lg1u16/R/win-library/3.5/jsonlite     | FALSE    | FALSE    | 2017-06-01 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| jtools       | jtools       | 1.1.1         | 1.1.1         | C:/Users/lg1u16/R/win-library/3.5/jtools       | C:/Users/lg1u16/R/win-library/3.5/jtools       | TRUE     | FALSE    | 2018-09-23 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| knitr        | knitr        | 1.20          | 1.20          | C:/Users/lg1u16/R/win-library/3.5/knitr        | C:/Users/lg1u16/R/win-library/3.5/knitr        | TRUE     | FALSE    | 2018-02-20 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| labeling     | labeling     | 0.3           | 0.3           | C:/Users/lg1u16/R/win-library/3.5/labeling     | C:/Users/lg1u16/R/win-library/3.5/labeling     | FALSE    | FALSE    | 2014-08-23 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| lattice      | lattice      | 0.20.35       | 0.20-35       | C:/Program Files/R/R-3.5.1/library/lattice     | C:/Program Files/R/R-3.5.1/library/lattice     | FALSE    | FALSE    | 2017-03-25 | CRAN (R 3.5.1)                         | NA    | C:/Program Files/R/R-3.5.1/library |
| lazyeval     | lazyeval     | 0.2.1         | 0.2.1         | C:/Users/lg1u16/R/win-library/3.5/lazyeval     | C:/Users/lg1u16/R/win-library/3.5/lazyeval     | FALSE    | FALSE    | 2017-10-29 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| lubridate    | lubridate    | 1.7.4         | 1.7.4         | C:/Users/lg1u16/R/win-library/3.5/lubridate    | C:/Users/lg1u16/R/win-library/3.5/lubridate    | FALSE    | FALSE    | 2018-04-11 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| magrittr     | magrittr     | 1.5           | 1.5           | C:/Users/lg1u16/R/win-library/3.5/magrittr     | C:/Users/lg1u16/R/win-library/3.5/magrittr     | FALSE    | FALSE    | 2014-11-22 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| MASS         | MASS         | 7.3.50        | 7.3-50        | C:/Program Files/R/R-3.5.1/library/MASS        | C:/Program Files/R/R-3.5.1/library/MASS        | TRUE     | FALSE    | 2018-04-30 | CRAN (R 3.5.1)                         | NA    | C:/Program Files/R/R-3.5.1/library |
| Matrix       | Matrix       | 1.2.14        | 1.2-14        | C:/Program Files/R/R-3.5.1/library/Matrix      | C:/Program Files/R/R-3.5.1/library/Matrix      | FALSE    | FALSE    | 2018-04-13 | CRAN (R 3.5.1)                         | NA    | C:/Program Files/R/R-3.5.1/library |
| memoise      | memoise      | 1.1.0         | 1.1.0         | C:/Users/lg1u16/R/win-library/3.5/memoise      | C:/Users/lg1u16/R/win-library/3.5/memoise      | FALSE    | FALSE    | 2017-04-21 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| modelr       | modelr       | 0.1.2         | 0.1.2         | C:/Users/lg1u16/R/win-library/3.5/modelr       | C:/Users/lg1u16/R/win-library/3.5/modelr       | FALSE    | FALSE    | 2018-05-11 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| MuMIn        | MuMIn        | 1.42.1        | 1.42.1        | C:/Users/lg1u16/R/win-library/3.5/MuMIn        | C:/Users/lg1u16/R/win-library/3.5/MuMIn        | TRUE     | FALSE    | 2018-07-23 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| munsell      | munsell      | 0.5.0         | 0.5.0         | C:/Users/lg1u16/R/win-library/3.5/munsell      | C:/Users/lg1u16/R/win-library/3.5/munsell      | FALSE    | FALSE    | 2018-06-12 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| nlme         | nlme         | 3.1.137       | 3.1-137       | C:/Program Files/R/R-3.5.1/library/nlme        | C:/Program Files/R/R-3.5.1/library/nlme        | TRUE     | FALSE    | 2018-04-07 | CRAN (R 3.5.1)                         | NA    | C:/Program Files/R/R-3.5.1/library |
| patchwork    | patchwork    | 0.0.1         | 0.0.1         | C:/Users/lg1u16/R/win-library/3.5/patchwork    | C:/Users/lg1u16/R/win-library/3.5/patchwork    | TRUE     | FALSE    | 2018-07-20 | Github (<thomasp85/patchwork@7fb35b1>) | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| pillar       | pillar       | 1.3.0         | 1.3.0         | C:/Users/lg1u16/R/win-library/3.5/pillar       | C:/Users/lg1u16/R/win-library/3.5/pillar       | FALSE    | FALSE    | 2018-07-14 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| pkgbuild     | pkgbuild     | 1.0.2         | 1.0.2         | C:/Users/lg1u16/R/win-library/3.5/pkgbuild     | C:/Users/lg1u16/R/win-library/3.5/pkgbuild     | FALSE    | FALSE    | 2018-10-16 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| pkgconfig    | pkgconfig    | 2.0.2         | 2.0.2         | C:/Users/lg1u16/R/win-library/3.5/pkgconfig    | C:/Users/lg1u16/R/win-library/3.5/pkgconfig    | FALSE    | FALSE    | 2018-08-16 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| pkgload      | pkgload      | 1.0.1         | 1.0.1         | C:/Users/lg1u16/R/win-library/3.5/pkgload      | C:/Users/lg1u16/R/win-library/3.5/pkgload      | FALSE    | FALSE    | 2018-10-11 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| plyr         | plyr         | 1.8.4         | 1.8.4         | C:/Users/lg1u16/R/win-library/3.5/plyr         | C:/Users/lg1u16/R/win-library/3.5/plyr         | FALSE    | FALSE    | 2016-06-08 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| prettyunits  | prettyunits  | 1.0.2         | 1.0.2         | C:/Users/lg1u16/R/win-library/3.5/prettyunits  | C:/Users/lg1u16/R/win-library/3.5/prettyunits  | FALSE    | FALSE    | 2015-07-13 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| processx     | processx     | 3.2.0         | 3.2.0         | C:/Users/lg1u16/R/win-library/3.5/processx     | C:/Users/lg1u16/R/win-library/3.5/processx     | FALSE    | FALSE    | 2018-08-16 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| ps           | ps           | 1.2.0         | 1.2.0         | C:/Users/lg1u16/R/win-library/3.5/ps           | C:/Users/lg1u16/R/win-library/3.5/ps           | FALSE    | FALSE    | 2018-10-16 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| pscl         | pscl         | 1.5.2         | 1.5.2         | C:/Users/lg1u16/R/win-library/3.5/pscl         | C:/Users/lg1u16/R/win-library/3.5/pscl         | TRUE     | FALSE    | 2017-10-10 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| purrr        | purrr        | 0.2.5         | 0.2.5         | C:/Users/lg1u16/R/win-library/3.5/purrr        | C:/Users/lg1u16/R/win-library/3.5/purrr        | TRUE     | FALSE    | 2018-05-29 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| R6           | R6           | 2.3.0         | 2.3.0         | C:/Users/lg1u16/R/win-library/3.5/R6           | C:/Users/lg1u16/R/win-library/3.5/R6           | FALSE    | FALSE    | 2018-10-04 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| raster       | raster       | 2.8.6         | 2.8-6         | C:/Users/lg1u16/R/win-library/3.5/raster       | C:/Users/lg1u16/R/win-library/3.5/raster       | TRUE     | FALSE    | 2018-11-07 | Github (<rspatial/raster@589de57>)     | NA    | C:/Users/lg1u16/R/win-library/3.5  |
| RColorBrewer | RColorBrewer | 1.1.2         | 1.1-2         | C:/Users/lg1u16/R/win-library/3.5/RColorBrewer | C:/Users/lg1u16/R/win-library/3.5/RColorBrewer | FALSE    | FALSE    | 2014-12-07 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| Rcpp         | Rcpp         | 0.12.19       | 0.12.19       | C:/Users/lg1u16/R/win-library/3.5/Rcpp         | C:/Users/lg1u16/R/win-library/3.5/Rcpp         | FALSE    | FALSE    | 2018-10-01 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| readr        | readr        | 1.1.1         | 1.1.1         | C:/Users/lg1u16/R/win-library/3.5/readr        | C:/Users/lg1u16/R/win-library/3.5/readr        | TRUE     | FALSE    | 2017-05-16 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| readxl       | readxl       | 1.1.0         | 1.1.0         | C:/Users/lg1u16/R/win-library/3.5/readxl       | C:/Users/lg1u16/R/win-library/3.5/readxl       | FALSE    | FALSE    | 2018-04-20 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| remotes      | remotes      | 2.0.1         | 2.0.1         | C:/Users/lg1u16/R/win-library/3.5/remotes      | C:/Users/lg1u16/R/win-library/3.5/remotes      | FALSE    | FALSE    | 2018-10-19 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| reshape      | reshape      | 0.8.8         | 0.8.8         | C:/Users/lg1u16/R/win-library/3.5/reshape      | C:/Users/lg1u16/R/win-library/3.5/reshape      | FALSE    | FALSE    | 2018-10-23 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| reshape2     | reshape2     | 1.4.3         | 1.4.3         | C:/Users/lg1u16/R/win-library/3.5/reshape2     | C:/Users/lg1u16/R/win-library/3.5/reshape2     | FALSE    | FALSE    | 2017-12-11 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| rgdal        | rgdal        | 1.3.6         | 1.3-6         | C:/Users/lg1u16/R/win-library/3.5/rgdal        | C:/Users/lg1u16/R/win-library/3.5/rgdal        | FALSE    | FALSE    | 2018-10-16 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| rlang        | rlang        | 0.3.0.1       | 0.3.0.1       | C:/Users/lg1u16/R/win-library/3.5/rlang        | C:/Users/lg1u16/R/win-library/3.5/rlang        | FALSE    | FALSE    | 2018-10-25 | CRAN (R 3.5.1)                         | NA    | C:/Users/lg1u16/R/win-library/3.5  |
| rmarkdown    | rmarkdown    | 1.10          | 1.10          | C:/Users/lg1u16/R/win-library/3.5/rmarkdown    | C:/Users/lg1u16/R/win-library/3.5/rmarkdown    | FALSE    | FALSE    | 2018-06-11 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| rprojroot    | rprojroot    | 1.3.2         | 1.3-2         | C:/Users/lg1u16/R/win-library/3.5/rprojroot    | C:/Users/lg1u16/R/win-library/3.5/rprojroot    | FALSE    | FALSE    | 2018-01-03 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| rstudioapi   | rstudioapi   | 0.8           | 0.8           | C:/Users/lg1u16/R/win-library/3.5/rstudioapi   | C:/Users/lg1u16/R/win-library/3.5/rstudioapi   | FALSE    | FALSE    | 2018-10-02 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| rvest        | rvest        | 0.3.2         | 0.3.2         | C:/Users/lg1u16/R/win-library/3.5/rvest        | C:/Users/lg1u16/R/win-library/3.5/rvest        | FALSE    | FALSE    | 2016-06-17 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| scales       | scales       | 1.0.0         | 1.0.0         | C:/Users/lg1u16/R/win-library/3.5/scales       | C:/Users/lg1u16/R/win-library/3.5/scales       | FALSE    | FALSE    | 2018-08-09 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| sessioninfo  | sessioninfo  | 1.1.0         | 1.1.0         | C:/Users/lg1u16/R/win-library/3.5/sessioninfo  | C:/Users/lg1u16/R/win-library/3.5/sessioninfo  | FALSE    | FALSE    | 2018-09-25 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| sf           | sf           | 0.7.1         | 0.7-1         | C:/Users/lg1u16/R/win-library/3.5/sf           | C:/Users/lg1u16/R/win-library/3.5/sf           | TRUE     | FALSE    | 2018-10-24 | CRAN (R 3.5.1)                         | NA    | C:/Users/lg1u16/R/win-library/3.5  |
| sp           | sp           | 1.3.1         | 1.3-1         | C:/Users/lg1u16/R/win-library/3.5/sp           | C:/Users/lg1u16/R/win-library/3.5/sp           | TRUE     | FALSE    | 2018-06-05 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| spData       | spData       | 0.2.9.4       | 0.2.9.4       | C:/Users/lg1u16/R/win-library/3.5/spData       | C:/Users/lg1u16/R/win-library/3.5/spData       | FALSE    | FALSE    | 2018-09-15 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| stringi      | stringi      | 1.2.4         | 1.2.4         | C:/Users/lg1u16/R/win-library/3.5/stringi      | C:/Users/lg1u16/R/win-library/3.5/stringi      | FALSE    | FALSE    | 2018-07-20 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| stringr      | stringr      | 1.3.1         | 1.3.1         | C:/Users/lg1u16/R/win-library/3.5/stringr      | C:/Users/lg1u16/R/win-library/3.5/stringr      | TRUE     | FALSE    | 2018-05-10 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| survival     | survival     | 2.42.6        | 2.42-6        | C:/Users/lg1u16/R/win-library/3.5/survival     | C:/Users/lg1u16/R/win-library/3.5/survival     | FALSE    | FALSE    | 2018-07-13 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| testthat     | testthat     | 2.0.1         | 2.0.1         | C:/Users/lg1u16/R/win-library/3.5/testthat     | C:/Users/lg1u16/R/win-library/3.5/testthat     | FALSE    | FALSE    | 2018-10-13 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| tibble       | tibble       | 1.4.2         | 1.4.2         | C:/Users/lg1u16/R/win-library/3.5/tibble       | C:/Users/lg1u16/R/win-library/3.5/tibble       | TRUE     | FALSE    | 2018-01-22 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| tidyr        | tidyr        | 0.8.2         | 0.8.2         | C:/Users/lg1u16/R/win-library/3.5/tidyr        | C:/Users/lg1u16/R/win-library/3.5/tidyr        | TRUE     | FALSE    | 2018-10-28 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| tidyselect   | tidyselect   | 0.2.5         | 0.2.5         | C:/Users/lg1u16/R/win-library/3.5/tidyselect   | C:/Users/lg1u16/R/win-library/3.5/tidyselect   | FALSE    | FALSE    | 2018-10-11 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| tidyverse    | tidyverse    | 1.2.1         | 1.2.1         | C:/Users/lg1u16/R/win-library/3.5/tidyverse    | C:/Users/lg1u16/R/win-library/3.5/tidyverse    | TRUE     | FALSE    | 2017-11-14 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| units        | units        | 0.6.1         | 0.6-1         | C:/Users/lg1u16/R/win-library/3.5/units        | C:/Users/lg1u16/R/win-library/3.5/units        | FALSE    | FALSE    | 2018-09-21 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| usethis      | usethis      | 1.4.0         | 1.4.0         | C:/Users/lg1u16/R/win-library/3.5/usethis      | C:/Users/lg1u16/R/win-library/3.5/usethis      | FALSE    | FALSE    | 2018-08-14 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| viridisLite  | viridisLite  | 0.3.0         | 0.3.0         | C:/Users/lg1u16/R/win-library/3.5/viridisLite  | C:/Users/lg1u16/R/win-library/3.5/viridisLite  | FALSE    | FALSE    | 2018-02-01 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| withr        | withr        | 2.1.2         | 2.1.2         | C:/Users/lg1u16/R/win-library/3.5/withr        | C:/Users/lg1u16/R/win-library/3.5/withr        | FALSE    | FALSE    | 2018-03-15 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| xml2         | xml2         | 1.2.0         | 1.2.0         | C:/Users/lg1u16/R/win-library/3.5/xml2         | C:/Users/lg1u16/R/win-library/3.5/xml2         | FALSE    | FALSE    | 2018-01-24 | CRAN (R 3.5.0)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |
| yaml         | yaml         | 2.2.0         | 2.2.0         | C:/Users/lg1u16/R/win-library/3.5/yaml         | C:/Users/lg1u16/R/win-library/3.5/yaml         | FALSE    | FALSE    | 2018-07-25 | CRAN (R 3.5.1)                         | TRUE  | C:/Users/lg1u16/R/win-library/3.5  |

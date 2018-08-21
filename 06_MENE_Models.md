MENE Models
================

``` r
library(MASS)
library(raster)
library(sf)
library(tidyverse)
library(GGally)
library(broom)
library(MuMIn)
library(knitr)

# resolutions
rlns <- c("1km", "5km", "10km", "25km", "50km", "100km")

# study extent (for MENE this is England)
study_ext_sf <- st_read("~/DATA/ADMINISTRATIVE/gb_shapefile/GBR_adm1.shp", quiet = TRUE) %>% 
  filter(NAME_1 == "England") %>% 
  st_transform(crs = 3035)

study_ext <- as(study_ext_sf, "Spatial")
```

Data collating
--------------

We're making one big dataframe of the response (MENE visits) and covariates (landcover diversity, topographical variation and population density) at each of the analysis resolutions.

``` r
get_df <- function(rln) {
  # list all files of specified resolution
  fnames <- list.files("data", pattern = paste0("_", rln), full.names = TRUE)
  
  # stack them
  dat <- stack(fnames)
  
  # crop by study extent
  dat <- crop(dat, study_ext)
  
  # get into a dataframe
  df <- as.data.frame(dat, xy = TRUE) %>% 
    rename_all(funs(str_replace_all(., paste0("_", rln), ""))) %>% 
    na.omit %>% 
    mutate(resolution = rln)
  
  return(df)
}

# dataframe with obs for all data
df <- map_dfr(rlns, get_df) %>% 
  mutate(resolution = factor(resolution, 
                             levels = c("1km", 
                                        "5km", 
                                        "10km", 
                                        "25km", 
                                        "50km", 
                                        "100km")))
```

Data exploration
----------------

First, let's just see whether there is sufficient variability and n in each resolution to fit models.

``` r
ggplot(df, aes(x = mene)) + 
  geom_histogram() + 
  facet_wrap(~resolution, scales = "free")
```

![](06_MENE_Models_files/figure-markdown_github/mene_distribution-1.png)

``` r
df %>% 
  group_by(resolution) %>% 
  summarise(count = n())
```

    ## # A tibble: 6 x 2
    ##   resolution count
    ##   <fct>      <int>
    ## 1 1km         5153
    ## 2 5km         2373
    ## 3 10km        1120
    ## 4 25km         254
    ## 5 50km          78
    ## 6 100km         24

Based on the distributions, we can fit models for all resolutions, but for the larger two (50km and 100km), n is very small.

We can also inspect the spatial coverage of the MENE data.

``` r
ggplot(df, aes(x = x, y = y)) + 
  geom_raster() + 
  coord_equal() + 
  facet_wrap(~resolution) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank())
```

![](06_MENE_Models_files/figure-markdown_github/mene_spatial-1.png)

Sparse coverage at 1km and 5km, but not too geographically biased. There are a lot of records around London, and to some extent within National Parks. This may cause some issues, but will see.

Finally, lets look at the relationships between variables:

``` r
p <- list()
for(r in rlns) {
  p[[r]] <- ggpairs(df %>% filter(resolution == r) %>% select(-resolution, -x, -y), progress = FALSE)
}

p
```

    ## $`1km`

![](06_MENE_Models_files/figure-markdown_github/var_relationships-1.png)

    ## 
    ## $`5km`

![](06_MENE_Models_files/figure-markdown_github/var_relationships-2.png)

    ## 
    ## $`10km`

![](06_MENE_Models_files/figure-markdown_github/var_relationships-3.png)

    ## 
    ## $`25km`

![](06_MENE_Models_files/figure-markdown_github/var_relationships-4.png)

    ## 
    ## $`50km`

![](06_MENE_Models_files/figure-markdown_github/var_relationships-5.png)

    ## 
    ## $`100km`

![](06_MENE_Models_files/figure-markdown_github/var_relationships-6.png)

Seems that it's population and proportion of non-built land covers that are the strongest predictor at all resolutions, and that this effect increases with spatial resolution. LC diversity is next, with similar effect. Note that population density and proportion of non-built land covers are highly correlated at all resolutions except for 50km and 100km.

Statistical Models
------------------

To start, we will fit poisson GLMs with all three covariates and their first-order interactions as the full model.

``` r
fit_mod <- function(dat) {
  glm(mene ~ clc_prop + clc_shei + dem_mean + pop,
      data = dat,
      family = "poisson")
}

scale_cols <- function(x) {
  scale_this <- function(y) as.vector(scale(y))
  #x <- mutate(x, pop = log10(pop + 1))
  mutate_at(x, .vars = vars(-mene), .funs = funs(scale_this))
}

# get the CA function
source("../upscaling/R/ca_glm.R")

# get all of the required stats
get_mod_stats <- function(mod) {
  coeff <- coef(mod)[-1] %>% 
    data.frame %>% 
    as_tibble(rownames = "variable") %>% 
    rename(coef = 2)
  
  conf <- confint(mod)[-1,] %>% 
    data.frame %>% 
    as_tibble(rownames = "variable") %>% 
    rename(lci = X2.5.., uci = X97.5..)

  ca <- calc_commonality(mod)$CCTotalbyVar %>% 
    as_tibble(rownames = "variable")
  
  out <- inner_join(coeff, conf) %>% 
    inner_join(ca)

  return(out)
}

mod <- df %>% 
  group_by(resolution) %>% 
  nest() %>% 
  mutate(data_scaled = map(data, scale_cols),
         mod = map(data_scaled, fit_mod),
         mod_glance = map(mod, glance),
         mod_stats = map(mod, get_mod_stats),
         mod_resid = map(mod, resid),
         mod_fitted = map(mod, fitted))

mod_D2 <- mod %>%
  select(resolution, mod_glance) %>% 
  unnest() %>% 
  mutate(D2 = (null.deviance - deviance)/null.deviance) %>% 
  select(resolution, D2)

mod_stats <- mod %>% 
  select(resolution, mod_stats) %>% 
  unnest() %>% 
  inner_join(mod_D2)
```

TODO: Testing residuals with DHARMa showed overdispersion with strong patterns in the residuals for CLC and POP. Would be worth fitting a quadratic term for each of these. Even with population logged, there are some low population cells (even at 100km level) with high residual values.

How well do they fit? NB We are using *D*<sup>2</sup>, which is equivalent to *R*<sup>2</sup> (Guisan and Zimmermann 2000). Note also that due to the differences in sample size, *D*<sup>2</sup> values are not directly comparable.

``` r
ggplot(mod_D2, aes(x = resolution, y = D2, group = 1)) + 
  geom_point() + 
  stat_summary(fun.y=sum, geom="line")
```

![](06_MENE_Models_files/figure-markdown_github/unnamed-chunk-1-1.png)

Explanatory power of the model increases greatly between 1km and 5km resolutions (with increasing explanatory power as resolution increases, as to be expected).

What about the relative importances?

``` r
ggplot(mod_stats, aes(x = resolution, y = Total/D2, fill = variable)) + 
  scale_fill_viridis_d() + 
  geom_bar(stat = "identity", position = "stack")
```

![](06_MENE_Models_files/figure-markdown_github/unnamed-chunk-2-1.png)

Population explains the most variation at all resolutions. Proportion of non-built land covers explains a large amount of variation at all but the coarsest resolutions (50km, 100km). The importance of land-cover diversity increases with spatial resolution (decreasing again at 100km, but this probably related to high correlation with proportion and small sample size). Finally, topography does not have a huge importance in the models, but mean elevation explains the most variation at 1km, this then decreases until little to no variance is explained at 100km. Variation in topography explains the least variance at all resolutions.

What about the relationships?

``` r
ggplot(mod_stats, aes(x = resolution, y = coef)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  geom_errorbar(aes(ymin = lci, ymax = uci)) + 
  facet_wrap(~variable, nrow = 1)
```

![](06_MENE_Models_files/figure-markdown_github/unnamed-chunk-3-1.png)

Some scale dependencies. Probably need to check they are not entirely statistical artefacts.

Bright & dark spots
-------------------

Following Frei et al. (2018) and Cinner et al. (2016), we calculate bright spots and dark spots. We define bright/dark spots as those where the observed value is +/- 1SD (of observed values) from the expected value.

``` r
np <- st_read("~/DATA/ADMINISTRATIVE/national_parks_england/National_Parks_England.shp",
              quiet = TRUE) %>%
  st_transform(st_crs(study_ext_sf))
# 
# cities <- st_read("~/DATA/ADMINISTRATIVE/uk_cities/Major_Towns_and_Cities_December_2015_Boundaries.shp",
#                   quiet = TRUE) %>% 
#   st_transform(st_crs(study_ext_sf)) %>% 
#   st_centroid %>% 
#   filter(tcity15nm == "London")

aonb <- st_read("~/DATA/ADMINISTRATIVE/aonb_england/Areas_of_Outstanding_Natural_Beauty_England.shp",
                quiet = TRUE) %>%
  st_transform(st_crs(study_ext_sf))

mod_fitted <- mod %>% 
  select(resolution, data, mod_fitted, mod_resid) %>% 
  unnest() %>% 
  group_by(resolution) %>% 
  mutate(diff = mod_fitted - mene,
         bd_spots = case_when(diff < -sd(mene) ~ "Dark",
                               diff > sd(mene) ~ "Bright",
                               TRUE ~ "Neutral"))

ggplot() + 
  geom_raster(data = mod_fitted, 
              aes(x = x, y = y, fill = bd_spots)) + 
  #geom_sf(data = np, colour = "black", fill = NA) + 
  #geom_sf(data = cities, colour = "black") + 
  #geom_sf(data = aonb, colour = "black", fill = NA) + 
  geom_sf(data = study_ext_sf, fill = NA) + 
  facet_wrap(~resolution) + 
  scale_fill_manual(values = c("orange", "blue", "grey")) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank())
```

![](06_MENE_Models_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
ggplot(mod_fitted, aes(x = mod_fitted, y = mene)) + 
  geom_point(aes(colour = bd_spots)) + 
  scale_colour_manual(values = c("orange", "blue", "grey")) + 
  geom_smooth(method = "lm") + 
  facet_wrap(~resolution, scales = "free")
```

![](06_MENE_Models_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
group_by(mod_fitted, resolution) %>% 
  summarise(Bright = sum(bd_spots == "Bright"),
            Dark = sum(bd_spots == "Dark"),
            Total = n(),
            `Bright %` = (Bright/Total)*100,
            `Dark %` = (Dark/Total)*100) %>% 
  kable
```

| resolution |  Bright|  Dark|  Total|   Bright %|    Dark %|
|:-----------|-------:|-----:|------:|----------:|---------:|
| 1km        |      37|   436|   5153|  0.7180283|  8.461091|
| 5km        |      69|   156|   2373|  2.9077118|  6.573957|
| 10km       |      27|    66|   1120|  2.4107143|  5.892857|
| 25km       |       4|     7|    254|  1.5748031|  2.755905|
| 50km       |       2|     1|     78|  2.5641026|  1.282051|
| 100km      |       0|     0|     24|  0.0000000|  0.000000|

References
----------

Cinner, Joshua E., Cindy Huchery, M. Aaron MacNeil, Nicholas A.J. Graham, Tim R. McClanahan, Joseph Maina, Eva Maire, et al. 2016. “Bright spots among the world’s coral reefs.” *Nature* 535 (7612). Nature Publishing Group: 416–19. doi:[10.1038/nature18607](https://doi.org/10.1038/nature18607).

Frei, Barbara, Delphine Renard, Matthew G. E. Mitchell, Verena Seufert, Rebecca Chaplin-Kramer, Jeanine M. Rhemtulla, and Elena M. Bennett. 2018. “Bright spots in agricultural landscapes: Identifying areas exceeding expectations for multifunctionality and biodiversity.” *Journal of Applied Ecology*, no. February: 1–13. doi:[10.1111/1365-2664.13191](https://doi.org/10.1111/1365-2664.13191).

Guisan, Antoine, and Niklaus E Zimmermann. 2000. “Predictive habitat distribution models in ecology.” *Ecological Modelling* 135 (2–3): 147–86. doi:[10.1016/s0304-3800(00)00354-9](https://doi.org/10.1016/s0304-3800(00)00354-9).

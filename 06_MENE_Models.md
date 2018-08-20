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
df <- map_dfr(rlns, get_df)
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
    ##   <chr>      <int>
    ## 1 100km         24
    ## 2 10km        1120
    ## 3 1km         5153
    ## 4 25km         254
    ## 5 50km          78
    ## 6 5km         2373

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
  p[[r]] <- ggpairs(df %>% filter(resolution == r) %>% select(-resolution, -x, -y))
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

Seems that it's population that is the strongest predictor at all resolutions, and that this effect increases with spatial resolution. LC diversity is second, with similar effect.

Statistical Models
------------------

To start, we will fit poisson GLMs with all three covariates and their first-order interactions as the full model.

``` r
fit_mod <- function(dat) {
  glm.nb(mene ~ clc * pop, 
      data = dat)
}

scale_cols <- function(x) {
  scale_this <- function(y) as.vector(scale(y))
  x <- mutate(x, pop = log10(pop + 1))
  mutate_at(x, .vars = vars(clc, dem, pop), .funs = funs(scale_this))
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
         resolution = factor(resolution, levels = c("1km", "5km", "10km", "25km", "50km", "100km")))

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

What about the relationships?

``` r
ggplot(mod_stats, aes(x = resolution, y = coef)) + 
  geom_point() + 
  geom_hline(yintercept = 0) + 
  geom_errorbar(aes(ymin = lci, ymax = uci)) + 
  facet_wrap(~variable, nrow = 1)
```

![](06_MENE_Models_files/figure-markdown_github/unnamed-chunk-3-1.png)

Bright/dark spots? Plot residuals &lt; 10th (dark) and &gt; 90th (bright) percentiles.

``` r
np <- st_read("~/DATA/ADMINISTRATIVE/national_parks_england/National_Parks_England.shp", 
              quiet = TRUE) %>% 
  st_transform(st_crs(study_ext_sf))

cities <- st_read("~/DATA/ADMINISTRATIVE/uk_cities/Major_Towns_and_Cities_December_2015_Boundaries.shp",
                  quiet = TRUE) %>% 
  st_transform(st_crs(study_ext_sf)) %>% 
  st_centroid %>% 
  filter(tcity15nm == "London")

aonb <- st_read("~/DATA/ADMINISTRATIVE/aonb_england/Areas_of_Outstanding_Natural_Beauty_England.shp",
                quiet = TRUE) %>% 
  st_transform(st_crs(study_ext_sf))

mod_resid <- mod %>% 
  select(resolution, data, mod_resid) %>% 
  unnest() %>% 
  group_by(resolution) %>% 
  mutate(resid_bin = case_when(mod_resid < quantile(mod_resid, 0.1) ~ "Dark",
                               mod_resid > quantile(mod_resid, 0.9) ~ "Bright",
                               TRUE ~ "Neutral"))

ggplot() + 
  geom_raster(data = mod_resid %>% filter(!resolution %in% c("1km", "5km")), 
              aes(x = x, y = y, fill = resid_bin)) + 
  geom_sf(data = np, colour = "black", fill = NA) + 
  geom_sf(data = cities, colour = "black") + 
  geom_sf(data = aonb, colour = "black", fill = NA) + 
  geom_sf(data = study_ext_sf, fill = NA) + 
  facet_wrap(~resolution) + 
  scale_fill_manual(values = c("orange", "blue", "grey")) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank())
```

![](06_MENE_Models_files/figure-markdown_github/unnamed-chunk-4-1.png)

We probably want to get rid of cells which are &lt; certain threshold land (dark spots in 50km and 100km maps likely down to this rather than actual dark spots).

Some of the bright spots are within, or near National Parks or AONBs (black polygons). Bright spots around London - which corresponds with the clustering of records we saw around London in the raw data. Does this suggest that the relationship with population is not given enough weighting?

References
----------

Guisan, Antoine, and Niklaus E Zimmermann. 2000. “Predictive habitat distribution models in ecology.” *Ecological Modelling* 135 (2–3): 147–86. doi:[10.1016/s0304-3800(00)00354-9](https://doi.org/10.1016/s0304-3800(00)00354-9).

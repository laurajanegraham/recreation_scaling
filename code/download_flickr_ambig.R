library(flickr)
library(tidyverse)

kw_landscape <- c("nature", "landscape", 
                  "cultural landscape", "cultural land", "hill",
                  "mountain", "valley", "basin", "highland",
                  "ridge", "cliff", "peak", "gorge",
                  "glacier", "beach", "shore", "coast",
                  "sea", "ocean", "wetland", "river", "dike",
                  "brook", "lake", "waterfall", "dune",
                  "swamp", "pond", "ditch", "channel",
                  "estuary", "creek", "forest", "tree",
                  "woods", "canopy", "grove", "hedgerow",
                  "bush", "meadow", "grassland", "pasture",
                  "countryside", "prairie", "maize",
                  "corn", "wheat", "oats", "livestock", "cattle",
                  "cow"," sheep", "orchard", "field",
                  "vineyard", "crops", "cropland", "grazing",
                  "heather", "heath", "heathland",
                  "park", "peat", "peatland", "peatbog",
                  "marsh", "marshes", "marshland",
                  "moor", "moors", "moorland", "shrubs",
                  "shrubland")

kw_ambig <- c("relax", "cruising", "relaxing", "beauty",
              "beautiful", "magnificence", "magnificent",
              "splendour", "brilliance", "brilliant",
              "inspiring", "inspired", "sublime",
              "gorgeous", "outstanding", "enjoying",
              "breathtaking", "enchanting")

kw_landscape <- str_sub(kw_landscape, " ", "")
kw_ambig <- str_sub(kw_ambig, " ", "")

jobid = as.integer(Sys.getenv("PBS_ARRAYID"))

kw <- kw_unambig[jobid]

p <- photo_by_tag(year_range = c(2009, 2017),
                  tags = kw,
                  woe_id = 24554868) %>% 
  mutate(photo_tags = str_split(photo_tags, pattern = " "),
         match = map_dbl(photo_tags, function(x) sum(x %in% kw_landscape))) %>% 
  filter(match > 0) %>% 
  select(-match)

save(p, file = paste0("data/flickr/", kw, ".Rda"))


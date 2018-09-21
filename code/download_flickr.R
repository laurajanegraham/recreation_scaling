library(flickr)

kw_unambig <- c("walk", "walking", "hike", "hiking",
                "camp", "camping", "recreation", "cycling",
                "horse riding", "fishing", "mountain biking", 
                "bike riding", "run", "running",
                "hunt", "hunting", "tourism",
                "climbing", "trekking", "mountaineering",
                "skiing", "sailing", "rowing", "jogging",
                "outdoor", "vista", "panorama",
                "scene", "scenic", "scenery", "view",
                "viewpoint", "heritage", "historic value")

jobid = as.integer(Sys.getenv("PBS_ARRAYID"))

kw <- kw_unambig[jobid]

p <- photosSearch(year_range = c(2009, 2017),
                  text = kw,
                  woe_id = 24554868)
save(p, file = paste0("data/flickr/", kw, ".Rda"))

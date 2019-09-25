library (spocc)
library (mapr)
library (dismo)
library(ENMeval)
library(raster)
library(rgdal)
library(rgbif)
library(maps)
library(dplyr)
library(ggplot2)

devtools::install_github('oshea-patrick/RSpatial')
library(RSpatial)

creosote_df = occ(query = "Larrea tridentata", 
                  from = "gbif",
                  limit = 22076,
                  has_coords = TRUE)
###############looking at predators of creosote
woodrat_df = occ(query = "Neotoma lepida", 
                   from = "gbif",
                   limit = 1000,
                   has_coords = TRUE)
ratoccdat = occ2df(woodrat_df)
ratloc=ratoccdat[,c('longitude', 'latitude')]

kangroo_df = occ(query = "Dipodomys", 
                 from = "gbif",
                 limit = 20000,
                 has_coords = TRUE)
kangoccdat = occ2df(kangroo_df)
kangoccdat = filter(kangoccdat, latitude < 40) %>% filter (latitude > 20) %>%
  filter(longitude > -125) %>% filter(longitude < -100 )
kangloc=kangoccdat[,c('longitude', 'latitude')]
kangloc = kangloc[-occ2thin,]

#######################################
map_leaflet(creosote_df)


occdat <- occ2df(creosote_df)
occdat_filtered = occdat
#Filter out outlier locations - all okay within first 1000 -- dplyer filter
occdat_filtered = filter(occdat, latitude < 40) %>% filter (latitude > 20) %>%
  filter(longitude > -125) %>% filter(longitude < -100 )



#sort(creosote_df2$date, decreasing = FALSE)


#clim = raster::stack(list.files('/usr/share/data/wc2.0/bio2.5/', pattern='.tif', full.names = TRUE))
#crop climate data down to north america



#plot data points on a climate map 

wc = getData('worldclim', var='bio', res = 5)



ext=extent(c(-125, -100, 20, 40))
wc2 = wc[[c('bio1', 'bio12', 'bio2', 'bio3')]] #annual factors (best so far)
wc3 = wc[[c('bio5', 'bio14', 'bio17', 'bio18')]] #Dryness and heat factors
wc4 = wc[[c('bio14', 'bio17', 'bio18')]] # extreme Dryness factors
wc5 = wc[[c('bio5', 'bio9', 'bio10')]] #extreme heat factors
wc6 = wc[[c('bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')]] #all rain factors
wc7 = wc[[c('bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11')]] #all temp factors
predictors = crop(wc2, ext)


loc=occdat[,c('longitude', 'latitude')]
extr = extract(predictors, loc)
loc = loc[!is.na(extr[,1]),]

#spatial thinning = more popular locations see more observations, even out the density of distribution
occ2thin = poThin(
  df = loc,
  spacing = 25,
  dimension = nrow(loc),
  lon = 'longitude',
  lat = 'latitude'
)

loc = loc[-occ2thin,]


eval = ENMevaluate(occ=as.data.frame(loc), env = predictors, method='block', parallel=FALSE, fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
best=which(eval@results$AICc == min(eval@results$AICc))

plot(eval@predictions[[best]])
points(as.data.frame(loc), pch=20, cex =0.1)

est.loc = extract(eval@predictions[[best]], as.data.frame(loc))
est.bg = extract(eval@predictions[[best]], eval@bg.pts)
ev = evaluate(est.loc, est.bg)
thr = threshold(ev)
plot(eval@predictions[[best]] > thr$equal_sens_spec, col = c('lightgrey', 'black'))
points(loc, pch = 20, cex= 1, col ='blue')
points(ratloc, pch = 20, cex= 1, col ='red')
points(kangloc, pch = 20, cex = 1, col = 'green')

#png('wc3.png')
#plot(eval@predictions[[best]] > thr$equal_sens_spec, col = c('lightgrey', 'black'))



###################
usmap = map_data("usa")
ggplot() + geom_polygon(data = usmap, 
                        aes(x=long, y = lat, group=group, color = "black"), fill= "white") + 
  coord_quickmap()  + geom_point(data = occdat_filtered, aes(x =longitude, y = latitude))


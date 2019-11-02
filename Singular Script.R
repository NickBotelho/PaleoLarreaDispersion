#Single Script
#Libraries####
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
library(RSpatial)
library(parallel)
library(maxnet)
devtools::install_github('oshea-patrick/RSpatial')


#GBIF pull####
creosote_df = occ(query = "Larrea tridentata", 
                  from = "gbif",
                  limit = 22076,
                  has_coords = TRUE)

occdat <- occ2df(creosote_df)
occdat_filtered = occdat

occdat_filtered = filter(occdat, latitude < 40) %>% filter (latitude > 20) %>%
  filter(longitude > -125) %>% filter(longitude < -100 )


diva_df = occ(query = "Larrea divaricata", 
              from = "gbif",
              limit = 1377,
              has_coords = TRUE)

divaoccdat = occ2df(diva_df)
divaoccdat = filter(divaoccdat, latitude < 15) %>% filter (latitude > -60) %>%
  filter(longitude > -90) %>% filter(longitude < -35 )
divaloc=divaoccdat[,c('longitude', 'latitude')]

#Load in Present and Paleo climate data (currently need hard coding)####
wc = getData('worldclim', var='bio', res = 5)
wc2 = wc[[c('bio1', 'bio12', 'bio2', 'bio3')]] #annual factors (best so far)

paleo_plots = raster::stack(list.files( pattern='.tif', full.names = TRUE)) #cclgmbi
paleo_plots2 = paleo_plots[[c('cclgmbi1','cclgmbi2','cclgmbi3','cclgmbi12')]] #raster using same variables

second_paleo_plots = raster::stack(list.files( pattern='.tif', full.names = TRUE)) #melgmbi
second_paleo_plots2 = second_paleo_plots[[c('melgmbi1','melgmbi2','melgmbi3','melgmbi12')]]

third_paleo_plots = raster::stack(list.files( pattern='.tif', full.names = TRUE))  #mrlgmbi
third_paleo_plots2 = third_paleo_plots[[c('mrlgmbi1','mrlgmbi2','mrlgmbi3','mrlgmbi12')]]

names(paleo_plots2) = c('bio1', 'bio2', 'bio3', 'bio12')
names(second_paleo_plots2) = c('bio1', 'bio2', 'bio3', 'bio12')
names(third_paleo_plots2) = c('bio1', 'bio2', 'bio3', 'bio12')

#Extents####
ext=extent(c(-125, -85, 10, 40)) #orginal ext southern US and mexico
paleo_ext = extent(c(-125, -20, -60, 40)) #Bounds to see from southern United states through Southern south america
paleo_central_america_ext = extent(c(-125, -25, -40, 30)) #closer look at central america

#spatial thinning = more popular locations see more observations, even out the density of distribution####
occ2thin = poThin(
  df = loc,
  spacing = 25,
  dimension = nrow(loc),
  lon = 'longitude',
  lat = 'latitude'
)
loc = loc[-occ2thin,]

#ENMeval for modern data####
#Tridentata
eval = ENMevaluate(occ=as.data.frame(loc), env = predictors, method='block', parallel=FALSE, fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
best=which(eval@results$AICc == min(eval@results$AICc))
#Divaracata
diva_eval = ENMevaluate(occ=as.data.frame(divaloc), env = predictors, method='block', parallel=FALSE, fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
diva_best=which(diva_eval@results$AICc == min(diva_eval@results$AICc))

#Calculate niche overlap####
trid_predict = predict(wc2, eval@models[[best]], type = "cloglog")
diva_predict = predict(wc2, diva_eval@models[[best]], type = "cloglog")

no = nicheOverlap(diva_predict, trid_predict, 
                  stat='D', mask=FALSE, checkNegatives=FALSE)
#determine significance
#take an equal random sample from both groups, build 2 models on the random points and calculate an overlap on that

noTest = function(tbl1, tbl2, clim){
  n1 = nrow(tbl1)
  n2 = nrow(tbl2)
  merge = rbind(tbl1, tbl2)
  s1 = sample(nrow(merge), size=n1)
  
  sub1 = merge[s1,]
  sub2 = merge[-s1,]
  
  #enmeval
  predictors = crop(wc2, paleo_ext)
  s1_eval = ENMevaluate(occ=as.data.frame(sub1), env = predictors, method='block', parallel=FALSE,
                        fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
  s2_eval = ENMevaluate(occ=as.data.frame(sub2), env = predictors, method='block', parallel=FALSE, 
                        fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
  
  
  s1_predict = predict(wc2, s1_eval@models[[best]], type = "cloglog")
  s2_predict = predict(wc2, s2_eval@models[[best]], type = "cloglog")
  
  no = nicheOverlap(s1_predict, s2_predict, 
                    stat='D', mask=FALSE, checkNegatives=FALSE) 
  return(no)
}
nicheOverlapStored = data.frame()
x1 = 1
for (x1 in 1:100){
  noTestResult = noTest(divaloc, loc, clim = clim)
  nicheOverlapStored = rbind(nicheOverlapStored, noTestResult)
}

#Plotting Into South America####
transfer = extent(c(-90, -35, -60, 15))
trans.wc2 = crop(wc2, paleo_ext)
trans.predict = predict(trans.wc2, eval@models[[best]], type = "cloglog")
plot(trans.predict)
points(divaloc, pch = 20, cex= 0.1, col ='black')

#model performance on L. divaricata
est.divaloc = extract(eval@predictions[[best]], as.data.frame(divaloc))
est.bg = extract(eval@predictions[[best]], eval@bg.pts)
ev.diva = evaluate(est.divaloc, est.bg)
thr = threshold(ev)


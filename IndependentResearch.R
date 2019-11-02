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


creosote_df = occ(query = "Larrea tridentata", 
                  from = "gbif",
                  limit = 22076,
                  has_coords = TRUE)

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

#extents
ext=extent(c(-125, -85, 10, 40)) #orginal ext southern US and mexico
paleo_ext = extent(c(-125, -20, -60, 40)) #Bounds to see from southern United states through Southern south america
paleo_central_america_ext = extent(c(-125, -25, -40, 30)) #closer look at central america
##
wc2 = wc[[c('bio1', 'bio12', 'bio2', 'bio3')]] #annual factors (best so far)
wc3 = wc[[c('bio5', 'bio14', 'bio17', 'bio18')]] #Dryness and heat factors
wc4 = wc[[c('bio14', 'bio17', 'bio18')]] # extreme Dryness factors
wc5 = wc[[c('bio5', 'bio9', 'bio10')]] #extreme heat factors
wc6 = wc[[c('bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')]] #all rain factors
wc7 = wc[[c('bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11')]] #all temp factors
predictors = crop(wc2, paleo_ext)


loc=occdat_filtered[,c('longitude', 'latitude')]
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

#Eval for the Tridentata
eval = ENMevaluate(occ=as.data.frame(loc), env = predictors, method='block', parallel=FALSE, fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
best=which(eval@results$AICc == min(eval@results$AICc))
######################################

########################################
plot(eval@predictions[[best]]) #green is high model score
points(as.data.frame(loc), pch=20, cex =0.1)

est.loc = extract(eval@predictions[[best]], as.data.frame(loc))
est.bg = extract(eval@predictions[[best]], eval@bg.pts)
ev = evaluate(est.loc, est.bg)
thr = threshold(ev)
plot(eval@predictions[[best]] > thr$equal_sens_spec, col = c('lightgrey', 'black'))
points(loc, pch = 20, cex= 0.1, col ='blue')


#plot(eval@predictions[[best]] > thr$equal_sens_spec, col = c('lightgrey', 'black'))



###################Plotting Into South America
diva_df = occ(query = "Larrea divaricata", 
              from = "gbif",
              limit = 1377,
              has_coords = TRUE)

divaoccdat = occ2df(diva_df)
divaoccdat = filter(divaoccdat, latitude < 15) %>% filter (latitude > -60) %>%
  filter(longitude > -90) %>% filter(longitude < -35 )
divaloc=divaoccdat[,c('longitude', 'latitude')]
#divaloc = divaloc[-occ2thin,]

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


### calculate niche overlap between L. tridentata model and L. divaricata
#ENMeval >  best model for divaricata
diva_eval = ENMevaluate(occ=as.data.frame(divaloc), env = predictors, method='block', parallel=FALSE, fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
diva_best=which(diva_eval@results$AICc == min(diva_eval@results$AICc))
#ENMeval > best model for tridentata
best=which(eval@results$AICc == min(eval@results$AICc))


trid_predict = predict(wc2, eval@models[[best]], type = "cloglog")
diva_predict = predict(wc2, diva_eval@models[[best]], type = "cloglog")

no = nicheOverlap(diva_predict, trid_predict, 
                  stat='D', mask=FALSE, checkNegatives=FALSE) 
#determine significance
#take an equal random sample from both groups, build 2 models on the random points and calculate an overlap on that
head(no)

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
#noTest(divaloc, loc, clim = clim)
#function to plot 70% diva, trid, and both into the past and return plots (all 3 of the past), do each diva, trid, both one at a time
#
temp = eval@models[[best]]

predictors = crop(wc2, paleo_ext)
evalRandomSample = function(tbl1, clim, predictors){
  plotSize = floor(0.3*(nrow(tbl1)))
  s1 = sample_n(tbl1, size=plotSize)
  
  s1_eval = ENMevaluate(occ=as.data.frame(s1), env = predictors, method='block', parallel=FALSE,
                        fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
  # tempBest=s1_eval@models[[which(s1_eval@results$AICc == min(s1_eval@results$AICc))]]
  
  best=which(s1_eval@results$AICc == min(s1_eval@results$AICc))
  
  #train maxent on s1
  # For picking model parameters on the complete set
  best_param = s1_eval@results[best, 1]
  best_arr = strsplit(as.character(best_param), "_")
  
  rm = best_arr[[1]][length(best_arr[[1]])]
  
  fc1 = best_arr[[1]][1:(length(best_arr[[1]]) - 1)]
  
  maxmatr = rbind(s1_eval@occ.pts, s1_eval@bg.pts)
  pres = c(rep(1, nrow(s1_eval@occ.pts)), rep(0, nrow(s1_eval@bg.pts)))
  maxmatr = cbind(maxmatr, pres)
  
  maxextr = extract(wc2, maxmatr[, c('LON', 'LAT')])
  best_mod = maxnet(
    p = maxmatr[, 'pres'],
    data = as.data.frame(maxextr),
    maxnet.formula(
      p = maxmatr[, 'pres'],
      data = as.data.frame(maxextr),
      classes = stringr::str_to_lower(fc1)
    ),
    regmult = as.numeric(rm)
  )
  return(best_mod)
  # return(tempBest)
}
paleo_plots2_crop = crop(paleo_plots2, paleo_ext)
second_paleo_plots2_crop = crop(second_paleo_plots2, paleo_ext)
third_paleo_plots2_crop = crop (third_paleo_plots2, paleo_ext)

plotRandomSample =function(raster1, tempBest, crop, clim, name){
  
  s1.predict1 = predict(crop, tempBest, type = "cloglog")
  names(s1.predict1) = paste(name, " Layer ")
  #raster1 = addLayer(raster1, s1.predict1)
  return(s1.predict1)
}
x1 = 1
diva_cclgmbi = plotRandomSample(diva_cclgmbi, tempBest, paleo_plots2_crop, clim = clim, "cclgmbi")
#Wrapper Function
sub_paleopredict <- function(x){
  tempBest = evalRandomSample(divaloc, clim = clim, predictors)
  diva_cclgmbi = plotRandomSample(diva_cclgmbi, tempBest, paleo_plots2_crop, clim = clim, "cclgmbi")
  diva_melgmbi = plotRandomSample(diva_melgmbi, tempBest, second_paleo_plots2_crop, clim = clim, "melgmbi")
  diva_mrlgmbi = plotRandomSample(diva_mrlgmbi, tempBest, third_paleo_plots2_crop, clim = clim, "mrlgmbi")
  ret = list(x, diva_cclgmbi, diva_melgmbi, diva_mrlgmbi)
  return(ret)
}

##Parallel
p_works2 = raster()
x1 = 1

nclus = 4 #Number of cores to allocate to this job:
cl = makeCluster(nclus, type ='FORK'); #Create cluster of size n
p_works2 = parLapply(cl, 1:16, sub_paleopredict) #Run parts across cluster and collect results
stopCluster(cl)

diva_cclgmbi = p_works2[[1]][[2]]
diva_melgmbi = p_works2[[1]][[3]]
diva_mrlgmbi = p_works2[[1]][[4]]

for(i in 2:length(p_works2)){
  diva_cclgmbi = stack(diva_cclgmbi, p_works2[[i]][[2]])
  diva_melgmbi = stack(diva_melgmbi, p_works2[[i]][[3]])
  diva_mrlgmbi = stack(diva_mrlgmbi, p_works2[[i]][[4]])
}
#Results of Parallel for diva##
diva_cclgmbi_mean = mean(diva_cclgmbi)
plot(diva_cclgmbi_mean, main ="diva_cclgmbi_mean")
diva_cclgmbi_max = max(diva_cclgmbi)
plot(diva_cclgmbi_max, main = "diva_cclgmbi_max")

plot(diva_cclgmbi_max - diva_cclgmbi_mean)
diva_melgmbi_mean = mean(diva_melgmbi)
plot(diva_melgmbi_mean)

diva_mrlgmbi_mean = mean(diva_mrlgmbi)
plot(diva_mrlgmbi_mean)





#Tridentata parallel#######################################################################
x1 = 1
p_works2 = raster()
trid_cclgmbi = raster()
trid_melgmbi = raster()
trid_mrlgmbi = raster()
trid_wc2 = raster()

sub_paleopredict <- function(x){
  tempBest = evalRandomSample(loc, clim = clim, predictors)
  trid_cclgmbi = plotRandomSample(trid_cclgmbi, tempBest, paleo_plots2_crop, clim = clim, "trid_cclgmbi")
  trid_melgmbi = plotRandomSample(trid_melgmbi, tempBest, second_paleo_plots2_crop, clim = clim, "trid_melgmbi")
  trid_mrlgmbi = plotRandomSample(trid_mrlgmbi, tempBest, third_paleo_plots2_crop, clim = clim, "trid_mrlgmbi")
  trid_wc2 = plotRandomSample(trid_wc2, tempBest, predictors, clim= clim, "trid_wc2")
  ret = list(x, trid_cclgmbi, trid_melgmbi, trid_mrlgmbi, trid_wc2)
  return(ret)
}

nclus = 4 #Number of cores to allocate to this job:
cl = makeCluster(nclus, type ='FORK'); #Create cluster of size n
p_works2 = parLapply(cl, 1:8, sub_paleopredict) #Run parts across cluster and collect results
stopCluster(cl)
  trid_cclgmbi = p_works2[[1]][[2]]
  trid_melgmbi = p_works2[[1]][[3]]
  trid_mrlgmbi = p_works2[[1]][[4]]
  trid_wc2 = p_works2[[1]][[5]]
  stopCluster(cl)
  
  for(i in 2:length(p_works2)){
    trid_cclgmbi = stack(trid_cclgmbi, p_works2[[i]][[2]])
    trid_melgmbi = stack(trid_melgmbi, p_works2[[i]][[3]])
    trid_mrlgmbi = stack(trid_mrlgmbi, p_works2[[i]][[4]])
    trid_wc2 = stack(trid_wc2, p_works2[[i]][[5]])
  }

#Tridentata Results##
trid_cclgmbi_mean = mean(trid_cclgmbi)
plot(trid_cclgmbi_mean)




#Combined Parallel##
combinedloc = rbind(loc, divaloc)
combined_cclgmbi = raster()
combined_melgmbi = raster()
combined_mrlgmbi = raster()
combined_wc2 = raster()

sub_paleopredict <- function(x){
  tempBest = evalRandomSample(loc, clim = clim, predictors)
  combined_cclgmbi = plotRandomSample(combined_cclgmbi, tempBest, paleo_plots2_crop, clim = clim, "combined_cclgmbi")
  combined_melgmbi = plotRandomSample(combined_melgmbi, tempBest, second_paleo_plots2_crop, clim = clim, "combined_melgmbi")
  combined_mrlgmbi = plotRandomSample(combined_mrlgmbi, tempBest, third_paleo_plots2_crop, clim = clim, "combined_mrlgmbi")
  combined_wc2 = plotRandomSample(combined_wc2, tempBest, predictors, clim = clim, "combined_wc2")
  ret = list(x, combined_cclgmbi, combined_melgmbi, combined_mrlgmbi, combined_wc2)
  return(ret)
}

p_works2 = raster()

  nclus = 4 #Number of cores to allocate to this job:
  cl = makeCluster(nclus, type ='FORK'); #Create cluster of size n
  p_works2 = parLapply(cl, 1:8, sub_paleopredict) #Run parts across cluster and collect results
  combined_cclgmbi = p_works2[[1]][[2]]
  combined_melgmbi = p_works2[[1]][[3]]
  combined_mrlgmbi = p_works2[[1]][[4]]
  combined_wc2 = p_works2[[1]][[5]]
  stopCluster(cl)

  for(i in 2:length(p_works2)){
    combined_cclgmbi = stack(combined_cclgmbi, p_works2[[i]][[2]])
    combined_melgmbi = stack(combined_melgmbi, p_works2[[i]][[3]])
    combined_mrlgmbi = stack(combined_mrlgmbi, p_works2[[i]][[4]])
    combined_wc2 = stack(combined_wc2, p_works2[[i]][[5]])
  }

#Combined Results##
#Pulling paleo raster from the last glacial maximum (cclgmbi)
#I got pulled these by setting the working directory directly into the files and just saving the list
#The numbers correspond to the same variable setups as wc from above
paleo_plots = raster::stack(list.files( pattern='.tif', full.names = TRUE)) #cclgmbi
paleo_plots2 = paleo_plots[[c('cclgmbi1','cclgmbi2','cclgmbi3','cclgmbi12')]] #raster using same variables
paleo_plots3 = paleo_plots[[c('cclgmbi5','cclgmbi14','cclgmbi17','cclgmbi18')]] #raster using same variables
paleo_plots4 = paleo_plots[[c('cclgmbi14','cclgmbi17','cclgmbi18')]] #extreme dryness BEST
paleo_plots5 = paleo_plots[[c('cclgmbi5','cclgmbi9','cclgmbi10')]] #extreme heat
paleo_plots6 = paleo_plots[[c('cclgmbi13','cclgmbi14','cclgmbi15','cclgmbi17','cclgmbi16','cclgmbi18','cclgmbi19')]] 
paleo_plots7 = paleo_plots[[c('cclgmbi5','cclgmbi6','cclgmbi7','cclgmbi8','cclgmbi9','cclgmbi10','cclgmbi11')]]
##The other paleo maps: melgmbi and mrlgmbi
second_paleo_plots = raster::stack(list.files( pattern='.tif', full.names = TRUE)) #melgmbi
second_paleo_plots2 = second_paleo_plots[[c('melgmbi1','melgmbi2','melgmbi3','melgmbi12')]]
second_paleo_plots4 = second_paleo_plots[[c('melgmbi14','melgmbi17','melgmbi18')]]

third_paleo_plots = raster::stack(list.files( pattern='.tif', full.names = TRUE))  #mrlgmbi
third_paleo_plots2 = third_paleo_plots[[c('mrlgmbi1','mrlgmbi2','mrlgmbi3','mrlgmbi12')]]
third_paleo_plots4 = third_paleo_plots[[c('mrlgmbi14','mrlgmbi17','mrlgmbi18')]]


paleo_ext = extent(c(-125, -20, -60, 40)) #Bounds to see from southern United states through Southern south america
paleo_central_america_ext = extent(c(-125, -25, -40, 30)) #closer look at central america

names(paleo_plots2) = c('bio1', 'bio2', 'bio3', 'bio12')
paleo.paleoplots2 = crop(paleo_plots2, paleo_ext)
paleo.predict = predict(paleo.paleoplots2, eval@models[[best]], type = "cloglog")
plot(paleo.predict)
###using the two new paleo maps
names(second_paleo_plots2) = c('bio1', 'bio2', 'bio3', 'bio12')
paleo.second_paleoplots2 = crop(second_paleo_plots2, paleo_ext)
paleo.predict2 = predict(paleo.second_paleoplots2, eval@models[[best]], type = "cloglog")
plot(paleo.predict2)

names(third_paleo_plots2) = c('bio1', 'bio2', 'bio3', 'bio12')
paleo.third_paleoplots2 = crop(third_paleo_plots2, paleo_ext)
paleo.predict3 = predict(paleo.third_paleoplots2, eval@models[[best]], type = "cloglog")
plot(paleo.predict3)





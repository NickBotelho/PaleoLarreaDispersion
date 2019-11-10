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
library(ggpubr)
library(viridis)
library(gridExtra)
devtools::install_github('oshea-patrick/RSpatial')

#Check if existing GBIF data exists####
file_list = list.files()
diva_exists = FALSE
trid_exists = FALSE
rastersExist = FALSE
climateDataExists = TRUE
for(i in 1:length(file_list)){
  if(file_list[i] == "diva_df.csv"){
    divaloc = read.csv2("diva_df.csv")
    diva_exists = TRUE
  }
  if (file_list[i] == "creosote_df.csv"){
    loc = read.csv2("creosote_df.csv")
    trid_exists = TRUE
  }
}

#GBIF pull####
if (trid_exists == FALSE){
creosote_df = occ(query = "Larrea tridentata", 
                  from = "gbif",
                  limit = 22076,
                  has_coords = TRUE)

occdat <- occ2df(creosote_df)
occdat_filtered = occdat
occdat_filtered = filter(occdat, latitude < 40) %>% filter (latitude > 20) %>%
  filter(longitude > -125) %>% filter(longitude < -100 )
loc=occdat_filtered[,c('longitude', 'latitude')]

#spatial thinning = more popular locations see more observations, even out the density of distribution####
occ2thin = poThin(
  df = loc,
  spacing = 25,
  dimension = nrow(loc),
  lon = 'longitude',
  lat = 'latitude'
)
loc = loc[-occ2thin,]


write.csv(loc,'creosote_df.csv')
}

if(diva_exists = FALSE){
diva_df = occ(query = "Larrea divaricata", 
              from = "gbif",
              limit = 1377,
              has_coords = TRUE)

divaoccdat = occ2df(diva_df)
divaoccdat = filter(divaoccdat, latitude < 15) %>% filter (latitude > -60) %>%
  filter(longitude > -90) %>% filter(longitude < -35 )
divaloc=divaoccdat[,c('longitude', 'latitude')]

write.csv(divaloc,'diva_df.csv')
}


if (climateDataExists == FALSE){
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
}
#Extents####
ext=extent(c(-125, -85, 10, 40)) #orginal ext southern US and mexico
paleo_ext = extent(c(-125, -20, -60, 40)) #Bounds to see from southern United states through Southern south america
paleo_central_america_ext = extent(c(-125, -25, -40, 30)) #closer look at central america

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

#Creating random 70% rasters####
  #Create the crops
predictors = crop(wc2, paleo_ext)
paleo_plots2_crop = crop(paleo_plots2, paleo_ext)
second_paleo_plots2_crop = crop(second_paleo_plots2, paleo_ext)
third_paleo_plots2_crop = crop (third_paleo_plots2, paleo_ext)
#Functions declared
evalRandomSample = function(tbl1, clim, predictors){
  plotSize = floor(0.6*(nrow(tbl1)))
  #plotPercentage = runif(1, min = 0.2, max = 0.95)
  #plotSize = floor( (plotPercentage)*(nrow(tbl1)))
  s1 = sample_n(tbl1, size=plotSize)
  
  s1_eval = ENMevaluate(occ=as.data.frame(s1), env = predictors, method='block', parallel=FALSE,
                        fc=c("L", "LQ"), RMvalues=seq(0.5, 2, 0.5), rasterPreds=T)
  #tempBest=s1_eval@models[[which(s1_eval@results$AICc == min(s1_eval@results$AICc))]]
  
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
  
}
plotRandomSample =function(raster1, tempBest, crop, clim, name){
  
  s1.predict1 = predict(crop, tempBest, type = "cloglog")
  names(s1.predict1) = paste(name, " Layer ")
  #raster1 = addLayer(raster1, s1.predict1)
  return(s1.predict1)
}

#Running the functions (if you need to)
for(i in 1:length(file_list)){
  if(file_list[i] == "diva_cclgmbi.gri"){
    rastersExist = TRUE
  }
}

if (rastersExist == FALSE){

for (i in 1:50){
  tempBest = evalRandomSample(divaloc, clim = clim, predictors)
  nextLayer = plotRandomSample(diva_cclgmbi, tempBest, paleo_plots2_crop, clim = clim, "diva_cclgmbi")
  diva_cclgmbi = stack(diva_cclgmbi, nextLayer)
  nextLayer = plotRandomSample(diva_melgmbi, tempBest, second_paleo_plots2_crop, clim = clim, "diva_melgmbi")
  diva_melgmbi = stack(diva_melgmbi, nextLayer)
  nextLayer = plotRandomSample(diva_cclgmbi, tempBest, third_paleo_plots2_crop, clim = clim, "diva_mrlgmbi")
  diva_mrlgmbi = stack(diva_mrlgmbi, nextLayer)
  nextLayer = plotRandomSample(diva_cclgmbi, tempBest, predictors, clim = clim, "diva_wc2")
  diva_wc2 = stack(diva_wc2, nextLayer)
  counter = counter + 1
}

for (i in 1:50){
  tempBest = evalRandomSample(loc, clim = clim, predictors)
  nextLayer = plotRandomSample(trid_cclgmbi, tempBest, paleo_plots2_crop, clim = clim, "trid_cclgmbi")
  trid_cclgmbi = stack(trid_cclgmbi, nextLayer)
  nextLayer = plotRandomSample(trid_melgmbi, tempBest, second_paleo_plots2_crop, clim = clim, "trid_melgmbi")
  trid_melgmbi = stack(trid_melgmbi, nextLayer)
  nextLayer = plotRandomSample(trid_cclgmbi, tempBest, third_paleo_plots2_crop, clim = clim, "trid_mrlgmbi")
  trid_mrlgmbi = stack(trid_mrlgmbi, nextLayer)
  nextLayer = plotRandomSample(trid_cclgmbi, tempBest, predictors, clim = clim, "trid_wc2")
  trid_wc2 = stack(trid_wc2, nextLayer)
}

for (i in 1:50){
  tempBest = evalRandomSample(combinedloc, clim = clim, predictors)
  nextLayer = plotRandomSample(combined_cclgmbi, tempBest, paleo_plots2_crop, clim = clim, "combined_cclgmbi")
  combined_cclgmbi = stack(combined_cclgmbi, nextLayer)
  nextLayer = plotRandomSample(combined_melgmbi, tempBest, second_paleo_plots2_crop, clim = clim, "combined_melgmbi")
  combined_melgmbi = stack(combined_melgmbi, nextLayer)
  nextLayer = plotRandomSample(combined_cclgmbi, tempBest, third_paleo_plots2_crop, clim = clim, "combined_mrlgmbi")
  combined_mrlgmbi = stack(combined_mrlgmbi, nextLayer)
  nextLayer = plotRandomSample(combined_cclgmbi, tempBest, predictors, clim = clim, "combined_wc2")
  combined_wc2 = stack(combined_wc2, nextLayer)
}
  #saving rasters
  writeRaster(diva_cclgmbi, "diva_cclgmbi")
  writeRaster(diva_melgmbi, "diva_melgmbi")
  writeRaster(diva_mrlgmbi, "diva_mrlgmbi")
  writeRaster(diva_wc2, "diva_wc2")
  
  writeRaster(trid_cclgmbi, "trid_cclgmbi")
  writeRaster(trid_melgmbi, "trid_melgmbi")
  writeRaster(trid_mrlgmbi, "trid_mrlgmbi")
  writeRaster(trid_wc2, "trid_wc2")
  
  writeRaster(combined_cclgmbi, "combined_cclgmbi")
  writeRaster(combined_melgmbi, "combined_melgmbi")
  writeRaster(combined_mrlgmbi, "combined_mrlgmbi")
  writeRaster(combined_wc2, "combined_wc2")
}

#Plot the rasters####
rasterToPlot = function(raster, color = "viridis"){
  df = as.data.frame(mean(raster), xy = TRUE)
  plot = ggplot() + 
    geom_raster(data = df,
                aes(x = x, y = y, fill = layer))
  #remove needless text
  plot = plot + theme(legend.position = 'none',  #remove legend
                      axis.title = element_blank(), #remove axis title
                      axis.ticks = element_blank(), #remove axis ticks
                      axis.text = element_blank()) #remove axis text
  #color scale
  plot = plot + scale_fill_viridis_c(option = color)
  return(plot)
}
g_legend <- function(a.gplot){ 
  plot = a.gplot + theme(legend.position='right') +theme(legend.key.height = unit(1, 'in')) + theme(legend.key.width = unit(1, 'in')) + 
    labs(fill = "Location Probabiltiy")
  tmp <- ggplot_gtable(ggplot_build(plot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]]
  
  return(legend)
} 

plot_diva_cclgmbi = rasterToPlot(diva_cclgmbi)
plot_diva_melgmbi = rasterToPlot(diva_melgmbi)
plot_diva_mrlgmbi = rasterToPlot(diva_mrlgmbi)
plot_diva_wc2 = rasterToPlot(diva_wc2)

plot_trid_cclgmbi = rasterToPlot(trid_cclgmbi)
plot_trid_melgmbi = rasterToPlot(trid_melgmbi)
plot_trid_mrlgmbi = rasterToPlot(trid_mrlgmbi)
plot_trid_wc2 = rasterToPlot(trid_wc2)

plot_combined_cclgmbi = rasterToPlot(combined_cclgmbi)
plot_combined_melgmbi =rasterToPlot(combined_melgmbi)
plot_combined_mrlgmbi =rasterToPlot(combined_mrlgmbi)
plot_combined_wc2 =rasterToPlot(combined_wc2)

plot_legend = g_legend(plot_diva_wc2) #create the legend alone

lay <- rbind(c(21,21,21,21,21,13),
             c(14,1,2,3,4,13),
             c(15,5,6,7,8,13),
             c(16,9,10,11,12,13),
             c(NA,17,18,19,20,13))
diva_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "L. Divaricata"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
trid_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "L. Tridentata"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
combined_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "Combined"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
wc2_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "wc2"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
cclgmbi_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "cclgmbi"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
melgmbi_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "melgmbi"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
mrlgmbi_label = ggplot() + geom_text(aes(x =0 , y = 0, label = "mrlgmbi"), size = 6) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())
main_title = ggplot() + geom_text(aes(x =0 , y = 0, label = "Divaricata and Tridentata predictions in paleo and modern climates"), size = 10) + theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank())


all_plots = grid.arrange(plot_diva_wc2, plot_diva_cclgmbi, plot_diva_melgmbi, plot_diva_mrlgmbi,
                         plot_trid_wc2, plot_trid_cclgmbi, plot_trid_melgmbi, plot_trid_mrlgmbi,
                         plot_combined_wc2, plot_combined_cclgmbi, plot_combined_melgmbi, plot_combined_mrlgmbi,
                         plot_legend,
                         diva_label, trid_label, combined_label, wc2_label, cclgmbi_label, melgmbi_label, mrlgmbi_label, main_title,
                         layout_matrix = lay)
all_plots

ggsave("all_plots.png", plot = all_plots, device = "png", width = 18, height = 9, unit = "in")


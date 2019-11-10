#plots
plotAnalyze = function(raster){
  name = substitute(raster)
  #rMax = max(raster)
  rMean = mean(raster)
  #plot(rMax, main = paste(name, " Max"))
  plot(rMean, main = paste(name, "Mean"))
  #plot(rMax - rMean, main = paste(name, "Max - Mean"))
  
}
rasterToPlot = function(raster){
  df = as.data.frame(mean(raster), xy = TRUE)
  plot = ggplot() + 
    geom_raster(data = df, 
                aes(x = x, y = y, fill = layer))
  plot = plot + theme(legend.position = 'none')
  return(plot)
  
}
plotAnalyze(diva_cclgmbi)
plotAnalyze(diva_melgmbi)
plotAnalyze(diva_mrlgmbi)
plotAnalyze(diva_wc2)

diva_wc2_df = as.data.frame(mean(diva_wc2), xy = TRUE)
plot_diva_wc2 = ggplot() + 
  geom_raster(data = diva_wc2_df, 
              aes(x = x, y = y, fill = layer))
#turn the rasters into df
#create a ggplot map
#create a multipannel ggplot thing

plot_diva_cclgmbi = rasterToPlot(diva_cclgmbi)
plot_diva_melgmbi = rasterToPlot(diva_melgmbi)
plot_diva_mrlgmbi = rasterToPlot(diva_mrlgmbi)
plot_diva = ggarrange(plot_diva_wc2, plot_diva_cclgmbi, plot_diva_melgmbi, plot_diva_mrlgmbi, ncol = 4, nrow = 1)
ggsave("plots_diva.png",plot = plot_diva, device = "png", width = 9, height = 9)

plotAnalyze(trid_cclgmbi)
plotAnalyze(trid_melgmbi)
plotAnalyze(trid_mrlgmbi)
plotAnalyze(trid_wc2)

plotAnalyze(combined_cclgmbi)
plotAnalyze(combined_melgmbi)
plotAnalyze(combined_mrlgmbi)
plotAnalyze(combined_wc2)

plot_trid_cclgmbi = rasterToPlot(trid_cclgmbi)
plot_trid_melgmbi = rasterToPlot(trid_melgmbi)
plot_trid_mrlgmbi = rasterToPlot(trid_mrlgmbi)
plot_trid_wc2 = rasterToPlot(trid_wc2)

plot_combined_cclgmbi = rasterToPlot(combined_cclgmbi)
plot_combined_melgmbi = rasterToPlot(combined_melgmbi)
plot_combined_mrlgmbi = rasterToPlot(combined_mrlgmbi)
plot_combined_wc2 = rasterToPlot(combined_wc2)

plot_trid = ggarrange(plot_trid_wc2, plot_trid_cclgmbi, plot_trid_melgmbi, plot_trid_mrlgmbi, ncol = 4, nrow = 1)
plot_combined = ggarrange(plot_combined_wc2, plot_combined_cclgmbi, plot_combined_melgmbi, plot_combined_mrlgmbi, ncol = 4, nrow = 1)
ggarrange(plot_diva, plot_trid, plot_combined, ncol = 4, nrow = 3)

#saving rasters
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

#library(ggpubr)
#gsa = ggarrange()
#ggsave()

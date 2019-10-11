#Extra code
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
map_leaflet(creosote_df) #leaflet
################
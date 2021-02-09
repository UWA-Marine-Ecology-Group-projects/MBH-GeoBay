library(measurements)
d<-read.csv("~/MBHdesignGB/Design3/GB_design3NPZ_HPZ_TransectsNotClus.csv")
d
head(d)

ds <-d

coordinates(ds) <- ~x+y

ds
bathy
proj4string(ds) <- proj4string(bathy)
ds2 <- spTransform(ds, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
ds2

writeOGR( ds2, dsn="~/MBHdesignGB/Design3", layer=paste( "hires_points", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
  
head(dsp2)

# change the degree symbol to a space
dsp2$y = gsub('.', ' ', dsp2$y)
dsp2$x = gsub('.', ' ', dsp2$x)

# convert from decimal minutes to decimal degrees
dsp2$lat = measurements::conv_unit(dsp2$y, from = 'dec_deg', to = 'deg_dec_min')
dsp2$long = measurements::conv_unit(dsp2$x, from = 'dec_deg', to = 'deg_dec_min')



lr <- read.csv("~/MBHdesignGB/Design4_notClusTrans_LatLon (1).csv")
head(lr)
# convert from decimal minutes to decimal degrees
lr$y.1= measurements::conv_unit(lr$y.1, from = 'dec_deg', to = 'deg_dec_min')
lr$x.1 = measurements::conv_unit(lr$x.1, from = 'dec_deg', to = 'deg_dec_min')

write.csv(lr, "coordinates_lores.csv")

lrs <- lr

bathy
b

coordinates(lrs) <- ~x+y
proj4string(lrs) <- proj4string(bathy)
lrs2 <- spTransform(lrs, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
lrs2

writeOGR( lrs2, dsn="~/MBHdesignGB/Design3", layer=paste( "lowres_points", Sys.Date(), sep="_"), driver="ESRI Shapefile", overwrite_layer=TRUE)

ldf <- as.data.frame(lrs2)
head(ldf)

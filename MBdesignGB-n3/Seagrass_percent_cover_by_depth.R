### Plot data based on depth, estimate SD of seagrass% cover by depth ---




library( raster)
library( rgdal)
library( sp)


rm( list=ls())


w.dir <- "~/MBHdesignGB"
s.dir <- "~/MBHdesignGB/SpatialData" # spatial data folder
o.dir <- "~/MBHdesignGB/outputs"
d.dir <- "~/MBHdesignGB/data"
p.dir <- "~/MBHdesignGB/plots"



#### AUV ####

# load bathy data ----
bathy <- raster(paste(s.dir, "GB_CMR_bathy.tif", sep='/'))
plot(bathy)
proj4string(bathy) # +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
bathy # 36395 cells

# load seagrass data ----
sg <- read.csv(paste(d.dir, "Auv_zoning.csv", sep='/'))
str(sg)

sgsp <- sg
coordinates(sgsp) <- ~coords.x1+coords.x2
proj4string(sgsp) <- proj4string(bathy)
points(sgsp, add=T)

## Get bathy for each AUV point ---
sgb <- raster::extract(bathy, sgsp, sp=T)
sgb
h <- hist(sgb$GB_CMR_bathy)
max(sgb$GB_CMR_bathy)
min(sgb$GB_CMR_bathy)
breaks <- c(-43, -33, -23, -13)
hist(sgb$GB_CMR_bathy, breaks = breaks)

## Plot by depth ranges: 13-23, 23-33, 33-43

sgd <- as.data.frame(sgb)
str(sgd)
names(sgd)
# remove unneeded columns --
sgd <- sgd[,c(2,3,8,15,23)]
head(sgd)


# create column of depth ranges --
# https://stackoverflow.com/questions/50988447/add-column-with-values-depending-on-another-column-to-a-dataframe
sgd$depthrange <- ifelse(sgd$GB_CMR_bathy > -23 & sgd$GB_CMR_bathy <= -13, '13-23m',
                 ifelse(sgd$GB_CMR_bathy > -33 & sgd$GB_CMR_bathy <= -23, '23-33m',
                        ifelse(sgd$GB_CMR_bathy <= -33, '33-43m', 'other depth')))
head(sgd)
str(sgd)
sgd$depthrange <- as.factor(sgd$depthrange)
levels(sgd$depthrange)

## calculate means, sd and se based on depth ranges --

smean <- stats::aggregate(Seagrasses ~ depthrange, data=sgd, FUN=mean)
smean
ssd <- stats::aggregate(Seagrasses ~ depthrange, data=sgd, FUN=sd)
ssd

## se function --
se <- function(x) sd(x) / sqrt(length(x)) # Create own function

sse <- stats::aggregate(Seagrasses ~ depthrange, data=sgd, FUN=se)
sse

## join them in one df --

datas <- cbind(smean, ssd$Seagrasses, sse$Seagrasses)
datas

datas
str(datas)
names <- c("Depth_range","mean", "sd", "se")
names(datas) <- names
names(datas)

#################

#### PLOTS ####

################

library(ggplot2)
library(ggthemes)
library(extrafont)
library(broman) # for colors: https://kbroman.files.wordpress.com/2014/05/crayons.png


## Set colors for plotting --
# Need 4 colours, one for each zone:
# for color list: https://kbroman.files.wordpress.com/2014/05/crayons.png
blue <- brocolors("crayons")["Turquoise Blue"] # "#77dde7"
green <- brocolors("crayons")["Inchworm"] # "#b2ec5d"
red <- brocolors("crayons")["Wild Watermelon"] # "#fc6c85" 
yellow <- brocolors("crayons")["Sunglow"] # "#ffcf48"


theme_set(theme_bw())
p<-ggplot(data=datas, aes(x=Depth_range, y=mean, fill = Depth_range)) +
  geom_bar(stat="identity", color = "black") +
  geom_errorbar(aes(ymax = mean-se, ymin = mean+se), width = 0.2, cex = 1) +
  geom_errorbar(aes(ymax = mean-sd, ymin = mean+sd), width = 0.2, color = "blue") +
  scale_fill_manual(values = c("#77dde7", "#fc6c85","#b2ec5d", "#ffcf48")) +
  labs(title = "AUV", y = "Seagrass mean % cover", x = "Depth range") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title.x = element_text(size = 12, face="bold"), axis.title.y = element_text(size = 12, face="bold"), 
        axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12, face="bold"),
        title = element_text(size = 14, face= "bold"))
p

ggsave(paste(p.dir, "Sg-depth-AUV.png", sep='/'), plot=p, device = "png", scale = 1, dpi =300 )


### Other depth ranges -----

sgd$depthrange2 <- ifelse(sgd$GB_CMR_bathy > -20 & sgd$GB_CMR_bathy <= -13, '13-20m',
                         ifelse(sgd$GB_CMR_bathy > -28 & sgd$GB_CMR_bathy <= -20, '20-28m',
                                       ifelse(sgd$GB_CMR_bathy <= -28, '28-43m', 'other depth')))
                                            
head(sgd)
str(sgd)
sgd$depthrange2 <- as.factor(sgd$depthrange2)
levels(sgd$depthrange2)

## calculate means, sd and se based on depth ranges --

smean <- stats::aggregate(Seagrasses ~ depthrange2, data=sgd, FUN=mean)
smean
ssd <- stats::aggregate(Seagrasses ~ depthrange2, data=sgd, FUN=sd)
ssd

## se function --
se <- function(x) sd(x) / sqrt(length(x)) # Create own function

sse <- stats::aggregate(Seagrasses ~ depthrange2, data=sgd, FUN=se)
sse

## join them in one df --

datas <- cbind(smean, ssd$Seagrasses, sse$Seagrasses)
datas

datas
str(datas)
names <- c("Depth_range","mean", "sd", "se")
names(datas) <- names
names(datas)

#################

#### PLOTS ####

################

library(ggplot2)
library(ggthemes)
library(extrafont)
library(broman) # for colors: https://kbroman.files.wordpress.com/2014/05/crayons.png


## Set colors for plotting --
# Need 4 colours, one for each zone:
# for color list: https://kbroman.files.wordpress.com/2014/05/crayons.png
blue <- brocolors("crayons")["Turquoise Blue"] # "#77dde7"
green <- brocolors("crayons")["Inchworm"] # "#b2ec5d"
red <- brocolors("crayons")["Wild Watermelon"] # "#fc6c85" 
yellow <- brocolors("crayons")["Sunglow"] # "#ffcf48"
yellowg <- brocolors("crayons")["Yellow Green"] #  "#c5e384"


theme_set(theme_bw())
p<-ggplot(data=datas, aes(x=Depth_range, y=mean, fill = Depth_range)) +
  geom_bar(stat="identity", color = "black") +
  geom_errorbar(aes(ymax = mean-se, ymin = mean+se), width = 0.2, cex = 1) +
  geom_errorbar(aes(ymax = mean-sd, ymin = mean+sd), width = 0.2, color = "blue") +
  scale_fill_manual(values = c("#77dde7", "#fc6c85","#b2ec5d", "#ffcf48", "#c5e384")) +
  labs(title = "AUV", y = "Seagrass mean % cover", x = "Depth range") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title.x = element_text(size = 12, face="bold"), axis.title.y = element_text(size = 12, face="bold"), 
        axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12, face="bold"),
        title = element_text(size = 14, face= "bold"))
p

ggsave(paste(p.dir, "Sg-depth-AUV2.png", sep='/'), plot=p, device = "png", scale = 1, dpi =300 )






#### TOWED VIDEO ####

rm( list=ls())


w.dir <- "~/MBHdesignGB"
s.dir <- "~/MBHdesignGB/SpatialData" # spatial data folder
o.dir <- "~/MBHdesignGB/outputs"
d.dir <- "~/MBHdesignGB/data"
p.dir <- "~/MBHdesignGB/plots"

# load bathy data ----
bathy <- raster(paste(s.dir, "GB_CMR_bathy.tif", sep='/'))
plot(bathy)
proj4string(bathy) # +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 
bathy # 36395 cells

# load seagrass data ----
sg <- read.csv(paste(d.dir, "TV_zoning.csv", sep='/'))
str(sg)

sgsp <- sg
coordinates(sgsp) <- ~coords.x1+coords.x2
proj4string(sgsp) <- proj4string(bathy)
points(sgsp, add=T) # this includes point outside the CMR

## Get bathy for each AUV point ---
sgb <- raster::extract(bathy, sgsp, sp=T)
sgb
h <- hist(sgb$GB_CMR_bathy)
max(sgb$GB_CMR_bathy) # many NAs
min(sgb$GB_CMR_bathy)
breaks <- c(-43, -33, -23, -13)
hist(sgb$GB_CMR_bathy, breaks = breaks)

## Plot by depth ranges: 13-23, 23-33, 33-43

sgd <- as.data.frame(sgb)
str(sgd)
names(sgd)
# remove unneeded columns --
sgd <- sgd[,c(2,6,13,21)]
head(sgd)
# remove NAs --
str(sgd) # 2962 obs.
sgd <- na.omit(sgd)
str(sgd) # 1129 obs.

hist(sgd$GB_CMR_bathy)
max(sgd$GB_CMR_bathy) # -12
min(sgd$GB_CMR_bathy) # -26

# create column of depth ranges --
# https://stackoverflow.com/questions/50988447/add-column-with-values-depending-on-another-column-to-a-dataframe
sgd$depthrange1 <- ifelse(sgd$GB_CMR_bathy > -20 & sgd$GB_CMR_bathy <= -12, '12-20m',
                         ifelse(sgd$GB_CMR_bathy <= -20, '20-26m', 'other depth'))
                               
head(sgd)
str(sgd)
sgd$depthrange <- as.factor(sgd$depthrange)
levels(sgd$depthrange)

## calculate means, sd and se based on depth ranges --

smean <- stats::aggregate(Seagrasses ~ depthrange, data=sgd, FUN=mean)
smean
ssd <- stats::aggregate(Seagrasses ~ depthrange, data=sgd, FUN=sd)
ssd

## se function --
se <- function(x) sd(x) / sqrt(length(x)) # Create own function

sse <- stats::aggregate(Seagrasses ~ depthrange, data=sgd, FUN=se)
sse

## join them in one df --

datas <- cbind(smean, ssd$Seagrasses, sse$Seagrasses)
datas

datas
str(datas)
names <- c("Depth_range","mean", "sd", "se")
names(datas) <- names
names(datas)

#################

#### PLOTS ####

################

library(ggplot2)
library(ggthemes)
library(extrafont)
library(broman) # for colors: https://kbroman.files.wordpress.com/2014/05/crayons.png


## Set colors for plotting --
# Need 4 colours, one for each zone:
# for color list: https://kbroman.files.wordpress.com/2014/05/crayons.png
blue <- brocolors("crayons")["Turquoise Blue"] # "#77dde7"
green <- brocolors("crayons")["Inchworm"] # "#b2ec5d"
red <- brocolors("crayons")["Wild Watermelon"] # "#fc6c85" 
yellow <- brocolors("crayons")["Sunglow"] # "#ffcf48"


theme_set(theme_bw())
p1<-ggplot(data=datas, aes(x=Depth_range, y=mean, fill = Depth_range)) +
  geom_bar(stat="identity", color = "black") +
  geom_errorbar(aes(ymax = mean-se, ymin = mean+se), width = 0.2, cex = 1) +
  geom_errorbar(aes(ymax = mean-sd, ymin = mean+sd), width = 0.2, color = "blue") +
  scale_fill_manual(values = c("#77dde7", "#fc6c85","#b2ec5d", "#ffcf48")) +
  labs(title = "TV", y = "Seagrass mean % cover", x = "Depth range") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title.x = element_text(size = 12, face="bold"), axis.title.y = element_text(size = 12, face="bold"), 
        axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12, face="bold"),
        title = element_text(size = 14, face= "bold"))
p1

ggsave(paste(p.dir, "Sg-depth-TV.png", sep='/'), plot=p1, device = "png", scale = 1, dpi =300 )

## Another depth range ----

# create column of depth ranges --
# https://stackoverflow.com/questions/50988447/add-column-with-values-depending-on-another-column-to-a-dataframe
sgd$depthrange2 <- ifelse(sgd$GB_CMR_bathy > -17 & sgd$GB_CMR_bathy <= -12, '12-17',
                          ifelse(sgd$GB_CMR_bathy > -21 & sgd$GB_CMR_bathy <= -17, '17-21m',
                          ifelse(sgd$GB_CMR_bathy <= -21, '21-26m', 'other depth')))

head(sgd)
str(sgd)
sgd$depthrange2 <- as.factor(sgd$depthrange2)
levels(sgd$depthrange2)

## calculate means, sd and se based on depth ranges --

smean <- stats::aggregate(Seagrasses ~ depthrange2, data=sgd, FUN=mean)
smean
ssd <- stats::aggregate(Seagrasses ~ depthrange2, data=sgd, FUN=sd)
ssd

## se function --
se <- function(x) sd(x) / sqrt(length(x)) # Create own function

sse <- stats::aggregate(Seagrasses ~ depthrange2, data=sgd, FUN=se)
sse

## join them in one df --

datas <- cbind(smean, ssd$Seagrasses, sse$Seagrasses)
datas

datas
str(datas)
names <- c("Depth_range","mean", "sd", "se")
names(datas) <- names
names(datas)

#################

#### PLOTS ####

################

library(ggplot2)
library(ggthemes)
library(extrafont)
library(broman) # for colors: https://kbroman.files.wordpress.com/2014/05/crayons.png


## Set colors for plotting --
# Need 4 colours, one for each zone:
# for color list: https://kbroman.files.wordpress.com/2014/05/crayons.png
blue <- brocolors("crayons")["Turquoise Blue"] # "#77dde7"
green <- brocolors("crayons")["Inchworm"] # "#b2ec5d"
red <- brocolors("crayons")["Wild Watermelon"] # "#fc6c85" 
yellow <- brocolors("crayons")["Sunglow"] # "#ffcf48"


theme_set(theme_bw())
p2<-ggplot(data=datas, aes(x=Depth_range, y=mean, fill = Depth_range)) +
  geom_bar(stat="identity", color = "black") +
  geom_errorbar(aes(ymax = mean-se, ymin = mean+se), width = 0.2, cex = 1) +
  geom_errorbar(aes(ymax = mean-sd, ymin = mean+sd), width = 0.2, color = "blue") +
  scale_fill_manual(values = c("#77dde7", "#fc6c85","#b2ec5d", "#ffcf48")) +
  labs(title = "TV", y = "Seagrass mean % cover", x = "Depth range") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
        axis.title.x = element_text(size = 12, face="bold"), axis.title.y = element_text(size = 12, face="bold"), 
        axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12, face="bold"),
        title = element_text(size = 14, face= "bold"))
p2

ggsave(paste(p.dir, "Sg-depth-TV2.png", sep='/'), plot=p2, device = "png", scale = 1, dpi =300 )







library(terra)
library(lubridate)  #for ym() and ceiling_date()
model <- c("EC-EARTH3")
scenario <- c("ssp585")
startym <- "2018-01"
endym <- "2020-12"
path <- "/home/james/Documents/data/"

alignRast <- function(input, target, outer, classes)
{
  i.crop2 <- crop(input, ext(target)+rep(outer,4)) #initially crop to larger extent than target (i.e. buffer)
  if(classes) { i.rs <- resample(i.crop2, target, method='near') }  #if categorical map
  else { i.rs <- resample(i.crop2, target) }
  i.crop <- crop(i.rs, ext(target))
  i.mask <- mask(i.crop, target)
  return(i.mask)
}

#target raster
munis.r <- rast('data/raster/socecon/G3MGsmunis_r_latlon.tif') #from Cerrado-LC-Input.Rmd

fpin <- paste0(path,model,"-pr-",scenario,".nc")
dat <- rast(fpin)
dat <- subset(dat,time(dat) >= ym(startym) & time(dat) < ceiling_date(ym(endym), unit="months"))

#tapp then align
start1 <- Sys.time()
datta <- tapp(dat, "yearmonths", sum, na.rm=T) #na.rm to handle 29 Feb missing data
datta <- terra::project(datta, "EPSG:4326")
datta <- alignRast(datta, munis.r, 2, TRUE)
end1 <- Sys.time() - start1
#Time difference of 10.94452 mins
#align then tapp
start2 <- Sys.time()
datat <- terra::project(dat, "EPSG:4326")
datat <- alignRast(datat, munis.r, 2, TRUE)
datat <- tapp(datat, "yearmonths", sum, na.rm=T) #na.rm to handle 29 Feb missing d
end2 <- Sys.time() - start2
#Time difference of 22.75285 mins
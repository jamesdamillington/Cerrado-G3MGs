#script to take rasters used with CRAFTY Brazil and align and output rasters for G3MGs study area

library(terra)

#align origin, resolution, extent and data cells (mask) of input raster to output raster 
#assumes crs of input and target are identical 
alignRast <- function(input, target, outer, classes)
{
  i.crop2 <- crop(input, ext(target)+rep(outer,4)) #initially crop to larger extent than target (i.e. buffer)
  if(classes) { i.rs <- resample(i.crop2, target, method='near') } 
  else { i.rs <- resample(i.crop2, target) }
  i.crop <- crop(i.rs, ext(target))
  i.mask <- mask(i.crop, target)
  return(i.mask)
}

#target raster
munis.r <- rast('data/raster/socecon/G3MGsmunis_r_latlon.tif') #from Cerrado-LC-Input.Rmd

#Transport capital
#input files created by https://github.com/jamesdamillington/CRAFTYInput/blob/master/PortAccessMap.r
rd2000 <- rast("data/raster/socecon/Transport/PortAccessCap2020_2000.asc")
rd2005 <- rast("data/raster/socecon/Transport/PortAccessCap2020_2005.asc")
rd2010 <- rast("data/raster/socecon/Transport/PortAccessCap2020_2010.asc")
rd2017 <- rast("data/raster/socecon/Transport/PortAccessCap2020_2017.asc")

rd2000 <- terra::project(rd2000, "EPSG:4326")
rd2005 <- terra::project(rd2005, "EPSG:4326")
rd2010 <- terra::project(rd2010, "EPSG:4326")
rd2017 <- terra::project(rd2017, "EPSG:4326")

cer2000 <- alignRast(rd2000, munis.r, 2, FALSE)
cer2005 <- alignRast(rd2005, munis.r, 2, FALSE)
cer2010 <- alignRast(rd2010, munis.r, 2, FALSE)
cer2017 <- alignRast(rd2017, munis.r, 2, FALSE)

writeRaster(cer2000, "data/raster/socecon/Transport/PortAccessCap2000_G3MGs.tif", "overwrite"=T)
writeRaster(cer2005, "data/raster/socecon/Transport/PortAccessCap2005_G3MGs.tif", "overwrite"=T)
writeRaster(cer2010, "data/raster/socecon/Transport/PortAccessCap2010_G3MGs.tif", "overwrite"=T)
writeRaster(cer2017, "data/raster/socecon/Transport/PortAccessCap2017_G3MGs.tif", "overwrite"=T)




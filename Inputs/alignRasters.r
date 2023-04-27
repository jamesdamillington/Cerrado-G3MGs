#script to take rasters used with CRAFTY Brazil and align and output rasters for G3MGs study area

library(terra)

#align origin, resolution, extent and data cells (mask) of input raster to output raster 
#assumes crs of input and target are identical 
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

#Protection Capital
#https://github.com/jamesdamillington/CRAFTYInput/blob/master/LandProtectionMap.r
protection <- rast("data/raster/socecon/LandProtection/All_ProtectionMap_025.asc")
protection <- terra::project(protection, "EPSG:4326")
protection <- alignRast(protection, munis.r, 2, TRUE)
writeRaster(protection, "data/raster/socecon/LandProtection/All_ProtectionMap_025_G3MGs.tif", "overwrite"=T)

#Land Value Capital
#input files created by https://github.com/jamesdamillington/CRAFTYInput/blob/master/LandValueMap.r
LV01 <- rast("data/raster/socecon/LandValue/LandPrice_Capital_08_2001.asc") 
LV17 <- rast("data/raster/socecon/LandValue/LandPrice_Capital_08_2017.asc") 

LV01 <- terra::project(LV01, "EPSG:4326")
LV17 <- terra::project(LV17, "EPSG:4326")

LV01 <- alignRast(LV01, munis.r, 2, TRUE)
LV17 <- alignRast(LV17, munis.r, 2, TRUE)

writeRaster(LV01, "data/raster/socecon/LandProtection/LandPrice_Capital_08_2001_G3MGs.tif", "overwrite"=T)
writeRaster(LV17, "data/raster/socecon/LandProtection/LandPrice_Capital_08_2017_G3MGs.tif", "overwrite"=T)

#Economic
econ <- rast("data/raster/socecon/Economic/AgriLocations2001.asc")
econ <- terra::project(econ, "EPSG:4326")
econ <- round(econ, 2)
econ <- alignRast(econ, munis.r, 2, TRUE)
writeRaster(econ, "data/raster/socecon/Economic/AgriLocations2001_G3MGs.tif", "overwrite"=T)

#Slope
slope <- rast("data/raster/physical/Slope/OAgri-slope_2018-08-16.asc")
slope <- terra::project(slope, "EPSG:4326")
hist(slope)  #check if should use categorical resample (answer: yes)
slope <- alignRast(slope, munis.r, 2, TRUE)
writeRaster(econ, "data/raster/physical/Slope/OAgri-slope_2018-08-16_G3MGs.tif", "overwrite"=T)

#HumanDev
hdev <- rast("data/raster/socecon/HumanDev/HumanCapital2001.asc")
hdev <- terra::project(hdev, "EPSG:4326")
hdev <- alignRast(hdev, munis.r, 2, TRUE)
writeRaster(hdev, "data/raster/socecon/HumanDev/HumanCapital2001_G3MGs.tif", "overwrite"=T)

#Soil
soil <- rast("data/raster/physical/soilT_2018-05-01.asc")
soil <- terra::project(soil, "EPSG:4326")
soil <- alignRast(soil, munis.r, 2, TRUE)
writeRaster(soil, "data/raster/physical/soilT_2018-05-01_G3MGs.tif", "overwrite"=T)
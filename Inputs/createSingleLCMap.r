#create maps of a single LC from observed LC maps
#useful for updating Other Agri and Other Capitals through simulation

#file copied from https://github.com/jamesdamillington/CRAFTYInput/blob/master/createSingleLCMap.r on 2022-07-09
#versioning tracked in https://github.com/jamesdamillington/Cerrado-G3MGs from 2022-07-09
#main edits are to 1) change package raster -> terra package 2) change to reclass all types and output multi-layer tif

rm(list=ls())
library(terra)

sim_yrs <- seq(2018, 2018, 1)   #single LC made for all these years
#sim_yrs <- c(2001, 2005, 2010)   #single LC made for all these years
reclass <- "reclass1"

#target names
tnames <- c("Nature", "OAgri", "Agri", "Other", "Pasture")  #order is target id

#if binary true output a binary map where 1 = target LC and 0 is not
#if binary false, output map where target LC value is maintained, all others set to 0
binary = T

obsLC <- rast(paste0("Inputs/data/raster/mapbiomas7/mapbiomas7-cerrado-G3MGs-2001-2021-1km-",reclass,".tif"))  #multi-layer

for(i in seq_along(sim_yrs)){
  
  outputLC  <- rast()
  
  for(target in seq_along(tnames)){
  
    #create df for subs below
    df <- data.frame(id=1:length(tnames), v=rep.int(0,length(tnames)))
    df[target,2] <- if(binary) 1 else target
    
    obsYr <- obsLC[paste0(sim_yrs[i])]  #paste0 to convert to character, terra uses partial matching of names (single []) 
    LC <- classify(obsYr, df)
    names(LC) <- tnames[target]

    outputLC <- c(outputLC, LC)  #warning message will be given because outputLC is initially empty
  }
  
  writeRaster(outputLC, paste0("Inputs/data/raster/mapbiomas7/mapbiomas7-cerrado-G3MGs-",sim_yrs[i],"-1km-",reclass,"-singleLCs.tif"), overwrite=T)
}

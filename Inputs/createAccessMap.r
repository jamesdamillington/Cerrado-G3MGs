#script to create 'access' capital maps

#requires maps output from createSingleLCMap.r 
#note: for singleLC_Nature_2001_PastureB_Disagg.asc this took approx 8 hours to complete (2.5 hours for Agri)

#file copied from https://github.com/jamesdamillington/CRAFTYInput/blob/master/accessMap.r on 2022-07-09
#versioning tracked in https://github.com/jamesdamillington/Cerrado-G3MGs from 2022-07-09
#main edits are to: 1. change raster -> terra package and 2. import multi-layer tifs (as createSingleLCMap.r is updated for terra)


library(raster)

#assumes focal value = 1
buf_width <- 5000  #width of buffer in m
BKGs <- list(0.05,0.05,0.0)  #background valeus
BUFs <- list(0.95,0.95,0.75) #buffer values
LCs <- list("Agri","OAgri","Nature") #the maps to work through
years <- c(2005, 2010)
suffix <- "_PastureB_Disagg.asc"


for(yr in years){
  
  for(i in seq_along(LCs)){
    
    print(yr)
    print(LCs[[i]])
    print(paste0("Start: ",Sys.time()))
    
    #create single LC map for the appropriate LC createSingleLCMap.r
    lc <- raster(paste0("Data/ObservedLCmaps/singleLC_",LCs[[i]],"_",yr,suffix))  #land cover from LandCoverMap.r (or ClassifyDisaggregateMap.r)
    
    crs(lc) <- latlong <- "+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs "
    
    #crop (for testing)
    #e <- extent(-60,-55,-15,-10)
    #lc <- crop(lc,e)
    #plot(lc)
    
    lc[lc == 0] <- NA  #set 0 to NA for buffer
    
    bf <- buffer(lc, width=buf_width, doEdge=T)  #buffer width of 5km
    
    s <- stack(lc, bf)
    
    s[is.na(s)] <- BKGs[[i]]
    #plot(s)
    
    setB <- function(a , b){
      ifelse(a != 1 & b == 1, BUFs[[i]], a)
    }
    
    out <- overlay(s, fun=setB)
    
    writeRaster(out, paste0("Data/",LCs[[i]],"Access_",yr,suffix), format = 'ascii', overwrite=T)
    
    print(paste0("End: ",Sys.time()))
  }
}
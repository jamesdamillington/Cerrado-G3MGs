#script to create 'access' capital maps

#requires maps output from createSingleLCMap.r 
#note: for singleLC_Nature_2001_PastureB_Disagg.asc this took approx 8 hours to complete (2.5 hours for Agri)

#file copied from https://github.com/jamesdamillington/CRAFTYInput/blob/master/accessMap.r on 2022-07-09
#versioning tracked in https://github.com/jamesdamillington/Cerrado-G3MGs from 2022-07-09
#main edits are to: 1. change raster -> terra package and 2. import multi-layer tifs (as createSingleLCMap.r is updated for terra)


library(terra)

#assumes focal value = 1
buf_width <- 1000  #width of buffer in m
#BKGs <- c(0.05,0.05,0.0)  #background valeusBKGs <- c(0.0,0.05,0.05,0.0,0.0)  #background valeus
#BUFs <- c(0.95,0.95,0.75) #buffer values
BUFs <- c(0.75,0.95,0.95,0.0,0.0) #buffer values
Cs <- c("Agri","OAgri","Nature") #the maps to work through
tnames <- c("Nature", "OAgri", "Agri", "Other", "Pasture")  #order is target id

targets <- c(1,2,3)

years <- c(2018)
reclass <- "reclass1"

yr <- 2018
#i<-2

for(yr in years){

  #singleLC (multi-layer) from createSingleLCMap.r
  #lcs <- rast(paste0("data/raster/mapbiomas6/mapbiomas6-cerrado-G3MGs-",yr,"-1km-",reclass,"-singleLCs.tif"))
  lc <- rast(paste0("mapbiomas7/mapbiomas7-cerrado-G3MGs-",yr,"-1km-",reclass,"-singleLCs.tif"))
  
  #reproject so that spatial unit is metres
  lcs <- project(lc, 'epsg:4087')
  
  for(i in targets){
    
    print(yr)
    print(tnames[i])
    print(paste0("Start: ",Sys.time()))
    
    lc <- lcs[[tnames[i]]]  #get this LC
    
    #crop (for testing)
    #e <- ext(-56,-55,-24,-23)
    #lc <- crop(lc,e)

    lc[lc == 0] <- NA  #set 0 to NA for buffer
    
    a <- global(lc, "notNA")  #needed for next check (any data in layer)
    
    if(a$notNA > 0) {  #if there is no cells of the class terra::buffer will crash
    
      bf <- buffer(lc, width=buf_width)  #"Unit is meter if x has a longitude/latitude CRS", from https://rspatial.github.io/terra/reference/buffer.html
      bf[!(bf)] <- NA
    } else { 
      bf <- lc     #if no cells 
    }
      
    s <- c(lc, bf)  #stack
    
    s[is.na(s)] <- BKGs[i]  #set cells that are neither target nor buffer to background value
    #plot(s)   
    
    #function sets target LC cells to 1, and buffer cells to buffer value 
    setB <- function(a , b){
      ifelse(a != 1 & b == 1, BUFs[i], a)
    }
    
    #apply the function
    out <- lapp(s, fun=setB)
    plot(out)

    writeRaster(out, paste0("data/raster/mapbiomas7/mapbiomas7-cerrado-G3MGs-",yr,"-1km-",reclass,"-",tnames[i],"Access.tif"), overwrite=T)
    
    print(paste0("End: ",Sys.time()))
  }
}

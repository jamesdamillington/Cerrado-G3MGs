#script to create moisture and GS capital maps from GCM output
#for climate data from ClimBRA Climate Change Dataset for Brazil https://doi.org/10.57760/sciencedb.02316
#based on moistureMap_future.r script for CRAFTY-Brazil that used  HadGEM2 data

library(terra)
library(tidyverse)
#library(ncdf4)
library(zoo)  #required for as.yearmon (masks terra::time)

#read munis.r as latlong
munis.r <- rast("data/raster/socecon/G3MGsmunis_r_latlon.tif") #map of municipality IDs to be simulated  #from alignRasters.r
munis.r

#soil is used to calculate plant available water
#read soil texture and set bucket size (see email from Daniel Victoria 2017-11-21)
#unzip(zipfile="Data/soilT_2018-05-01.zip",exdir="Data")  #unzip
soil<-raster("data/raster/physical/soilT_2018-05-01_G3MGs.tif")  
soil

#PAW is Plant Available Water
PAW<-soil
PAW[PAW==1]<-0
PAW[PAW==2]<-0
PAW[PAW==3]<-75
PAW[PAW==4]<-55
PAW[PAW==5]<-35
#plot(PAW)



model <- c("EC-EARTH3")
scenario <- c("ssp585")
startym <- "2015-01"
endym <- "2100-12"
path <- "/home/james/Documents/data/"


#sets time & names after reading monthly summarised ClimBRA data from ClimBRA-summar.r 
setClimBRAmonth <- function(crast, startlab, endlab){
  tc <- crast
  se <- seq(ym(startlab),ym(endlab),by="months")
  time(tc, tstep="yearmonths") <- se
  names(tc) <- mapply(paste, year(se), month(se, label=T))
  return(tc)
}



pr <- rast(paste0(path,model,"-pr-",scenario,"-totalmonth-",startym,"-to-",endym,".nc"))
pr <- setClimBRAmonth(pr, startym, endym) 

tasmin <- rast(paste0(path,model,"-tasmin-",scenario,"-meanmonth-",startym,"-to-",endym,".nc"))
tasmin <- setClimBRAmonth(tasmin, startym, endym) 

tasmax <- rast(paste0(path,model,"-tasmax-",scenario,"-meanmonth-",startym,"-to-",endym,".nc"))
tasmax <- setClimBRAmonth(tasmax, startym, endym) 




calcMoistureMaps <- function(munis.r, PAW, year, BRA.e, hemi, season, GS, RCP)
{
  #generate timeRange and filenames for this year
  
  
  #month labels for plots
  monlab <- c(paste0("Jan",year),paste0("Feb",year),paste0("Mar",year),paste0("Apr",year),paste0("May",year),paste0("Jun",year),paste0("Jul",year),paste0("Aug",year),paste0("Sep",year),paste0("Oct",year),paste0("Nov",year),paste0("Dec",year))
  if(hemi == "S")  monlab <- c(paste0("Jul",year-1),paste0("Aug",year-1),paste0("Sep",year-1),paste0("Oct",year-1),paste0("Nov",year-1),paste0("Dec",year-1),paste0("Jan",year),paste0("Feb",year),paste0("Mar",year),paste0("Apr",year),paste0("May",year),paste0("Jun",year))
  
  #create season filelabel
  season_label <- paste0(season, collapse="")
  if(hemi == "S" && season_label == "JulAugSepOctNovDecJanFebMarAprMayJun") season_label <- "All"
  if(hemi == "N" && season_label == "JanFebMarAprMayJunJulAugSepOctNovDec") season_label <- "All"
  
  
  #read climate files for this year 
  #N hemi is 1:12
  #S hemi 7:18
  #read climate files
  
  #northern hemisphere
  if(hemi == "N") {
    pre <- subset(pr,time(pr) >= year & time(pr) < year+1) 
    tmn <- subset(tasmin,time(tasmin) >= year & time(tasmin) < year+1) 
    tmx <- subset(tasmax,time(tasmax) >= year & time(tasmax) < year+1) 
  }
  
  #southern hemisphere
  if(hemi == "S") {
    pre <- subset(pr, time(pr) >= year+0.5 & time(pr) < year+1.5)
    tmn <- subset(tasmin, time(tasmin) >= year+0.5 & time(tasmin) < year+1.5)
    tmx <- subset(tasmax, time(tasmax) >= year+0.5 & time(tasmax) < year+1.5)
  }
  
  
  #project and crop bricks to extent we want for Cerrado
  
  #calculate average temperature by month (brick)
  avtemp.b <- 0.36*(3*tmx.b-tmn.b)

  #annual temp raster layer
  Ta<-mean(avtemp.b)
  
  #set params needed to calculate P and PET
  Days <- list(31,28,31,30,31,30,31,31,30,31,30,31) #list of days in the month 
  if(hemi == "S") {
    Days <- list(31,31,30,31,30,31,31,28,31,30,31,30)
  }
  
  #convert kg m-2 s-1 to mm/day
  pre.b <- pre.b * 86400 #s in day
  monthProd <- function(P, D) { return(P * D) }  #multiply by the number of days in the month
  pre.b <- 
    map2(as.list(pre.b), Days, monthProd) %>% 
    stack()  #remember to re-stack the list after function
  
  
  #total pptn raster layer
  Pa<-sum(pre.b)
  
  Idex <- 12*((0.2*Ta)^1.514)#1st regional thermal index
  Adex <- 0.49239+1.7912*10^-2*Idex-7.71*10^-5*Idex^2+6.75*10^-7*Idex^3#2nd regional thermal index
  Ndex <- 12##NEED REAL VALUE
  
  
  #initialize PET with mean temperatures
  PET.b <- avtemp.b 
  
  #function to calculate Potential Evapotranspiration (PET)
  calcPET <- function(PET, D, I, a, N)
  {
    #if mean temperature <= 0
    PET[PET<=0]<-0
    
    #if mean temperature > 26.5
    PET[PET>26.5]<-(-415.85+32.24*PET[PET>26.5]-0.43*(PET[PET>26.5]^2))*(N/12)*(D/30)
    
    #else
    PET[PET>0&PET<=26.5]<-0.0444*((10*(PET[PET>0&PET<=26.5]/I[PET[]>0&PET[]<=26.5]))^a[PET[]>0&PET[]<=26.5])*N*D
    
    return(PET)
  }
  
  
  #map2 to loop over raster brick (layer per month) and Days list (from purrr, see http://r4ds.had.co.nz/iteration.html)
  #brick needs to passed as a list (see https://geocompr.robinlovelace.net/location.html)
  PET.b <- 
    map2(as.list(PET.b), Days, calcPET, I = Idex, a = Adex, N = Ndex) %>% 
    stack()  #remember to re-stack the list after function
  
  names(PET.b) <- monlab
  
  #see Victoria et al. 2007 DOI: 10.1175/EI198.1 Table 2 for equations
  #initialise water storage variables 
  
  Stoi <- PAW  #Stoii is month i-1 storage
  Stoii <- Stoi #Stoi is month i storage (in first month use same values)
  
  allmeanStoi <- vector("double", 12)  #vector to hold meanStoi for each month
  
  #for creating empty rasters and bricks
  nullRaster <- munis.r
  nullRaster[!is.na(nullRaster)] <- 0  #set anywhere that is not 'no data' in munis.r to 0
  
  DEF.b <- stack(replicate(12, nullRaster)) #empty brick to save all month's DEF
  ET.b <- stack(replicate(12, nullRaster))  #empty brick to save all month's ET
  
  DEF <- nullRaster #empty layer for temp useage in loop
  ET <- nullRaster #empty layer for temp useage in loop
  
  #par(mfrow=c(1,1))
  
  #see loopProofs (need to use loop, cannot use map)
  for(i in 1:12)
  {
    #hold current values of Stoi to set Stoii for next month below (this is why we can't use map)
    tempStoi <- Stoi
    
    P <- pre.b[[i]]    #get this month's precipitation (for clarity in equations below)
    PET <- PET.b[[i]]  #get this month's PET (for clarity in equations below)
    
    #if pptn < PET set storage
    Stoi[P<PET] <- Stoii[P<PET] * exp(P[P<PET] - PET[P<PET]/PAW[P<PET])
    
    #if pptn >= PET set storage
    Stoi[P>=PET] <- Stoii[P>=PET] + (P[P>=PET] - PET[P>=PET])
    
    #update Stoii ready for next month
    Stoii<-tempStoi
    
    #where Sto > PAW
    Stoi[Stoi[]>PAW[]] <- PAW[Stoi[]>PAW[]]
    
    #save mean Stoi value for this month
    allmeanStoi[i] <- cellStats(Stoi, "mean")
    
    #calculate delta storage
    trSto <- Stoi - Stoii
    
    #reset ET for this loop
    ET <- nullRaster
    
    #where pptn < PET
    ET[P<PET] <- P[P<PET] - trSto[P<PET]
    
    #where P >= PET
    ET[P>=PET] <- PET[P>=PET]
    
    #reset DEF for this loop 
    DEF <- nullRaster
    
    #where pptn < PET
    DEF[P<PET] <- PET[P<PET] - ET[P<PET]
    
    #where P >= PET
    DEF[P>=PET]<-0
    
    #copy DEF to DEF brick
    DEF.b[[i]] <- DEF
    ET.b[[i]] <- ET
    
  }
  
  names(DEF.b) <- monlab    #apply month-year names to the layes (for plotting and access below)
  names(ET.b) <- monlab
  
  si_yr <- lapply(season, paste0, year)  #add the year to month names to match monlab format
  season_indices <- ifelse(monlab %in% si_yr,1,2)  #create index of months to use in stackApply below (1 is in season, 2 is not)
  
  #calculate Dryness Index
  #avDEF<-mean(DEF.b)#mean annual DEF
  avDEF <- stackApply(DEF.b, season_indices, mean)  #mean DEF (for specified months), #creates a stack of two layers (season months and non-season months)
  #avPET<-mean(PET.b)#mean annual PET
  avPET <- stackApply(PET.b, season_indices, mean)  #mean PET (for specified months), #creates a stack of two layers (season months and non-season months)
  
  avDi <- (100*avDEF) / avPET  #creates a stack of two layers (season months and non-season months)
  Di <- (100*DEF.b) / PET.b
  
  
  #pptn and temp by season if needed
  sumPptn <- stackApply(pre.b, season_indices, sum)
  avTemp <- stackApply(avtemp.b, season_indices, mean)
  
  
  #Number of months with water deficit - helper function
  countWD <- function(vect, na.rm=T) { return(sum(vect > 5)) }
  
  #Number of months in the season - helper function
  countMonths <- function(vect, na.rm=T) { return(length(vect)) }
  
  #calculate various ways of defining water deficif months
  DEFmonths <- stackApply(DEF.b,season_indices, countWD)  #creates a stack of two layers (season months and non-season months)
  allmonths <- stackApply(DEF.b,season_indices, countMonths) #creates a stack of two layers (season months and non-season months)
  DEFmonths_prop <- DEFmonths / allmonths   #proportion,  #creates a stack of two layers (season months and non-season months)
  
  #Stoidiffc <- allmeanStoi[12] - allmeanStoi[1]  #this does not seem to be used elsewhere...
  #Stoidiffc
  
  #write data to files
  if(writeClimRast)
  {
    writeRaster(avDEF[["index_1"]], paste0(outputDir,"/",className,"/MeanDEF_",season_label,"_","rcp",RCP,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(avPET[["index_1"]], paste0(outputDir,"/",className,"/MeanPET_",season_label,"_","rcp",RCP,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(avTemp[["index_1"]], paste0(outputDir,"/",className,"/MeanTemp_",season_label,"_","rcp",RCP,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(sumPptn[["index_1"]], paste0(outputDir,"/",className,"/SumPrecip_",season_label,"_","rcp",RCP,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(avDi[["index_1"]], paste0(outputDir,"/",className,"/MeanDI_",season_label,"_","rcp",RCP,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(DEFmonths[["index_1"]], paste0(outputDir,"/",className,"/CountDEFmonths_",season_label,"_","rcp",RCP,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
    writeRaster(DEFmonths_prop[["index_1"]], paste0(outputDir,"/",className,"/PropDEFmonths_",season_label,"_","rcp",RCP,"_",year,hemi,".asc"), format = 'ascii', overwrite=T)
  }
  
  #write pdfs
  if(writeClimPdf)
  {
    
    pdf(paste0(outputDir,"/",className,"/DEF_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(DEF.b, ext = BRA.e)
    dev.off()
    
    pdf(paste0(outputDir,"/",className,"/PET_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(PET.b, ext = BRA.e)
    dev.off()
    
    pdf(paste0(outputDir,"/",className,"/ET_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(ET.b, ext = BRA.e)
    dev.off()
    
    pdf(paste0(outputDir,"/",className,"/PPTNmon_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(pre.b, ext = BRA.e)
    dev.off()
    
    pdf(paste0(outputDir,"/",className,"/SumPrecip_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(sumPptn[["index_1"]], ext = BRA.e)
    dev.off()
    
    
    #pptnname <- ""
    #if(GS) { pptnname <- paste0(outputDir,"/",className,"/PPTN_GS","rcp",RCP,"_",year,hemi,".pdf") }
    #if(!GS) { pptnname <- paste0(outputDir,"/",className,"/PPTN_Mois","rcp",RCP,"_",year,hemi,".pdf") }
    
    #pdf(pptnname)
    #plot(sumPptn, ext = BRA.e)
    #dev.off()
    
    pdf(paste0(outputDir,"/",className,"/MeanTemp_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(avTemp[["index_1"]], ext = BRA.e)
    dev.off()
    
    #tempname <- ""
    #if(GS) { tempname <- paste0(outputDir,"/",className,"/TEMP_GS","rcp",RCP,"_",year,hemi,".pdf") }
    #if(!GS) { tempname <- paste0(outputDir,"/",className,"/TEMP_Mois","rcp",RCP,"_",year,hemi,".pdf") }
    
    #pdf(tempname)
    #plot(avTemp, ext = BRA.e)
    #dev.off()
    
    pdf(paste0(outputDir,"/",className,"/DI_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(Di, ext = BRA.e)
    dev.off()
    
    pdf(paste0(outputDir,"/",className,"/ClimateVariables_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    #pdf(paste0(outputDir,"/",className,"/meanDEF_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(avDEF[["index_1"]], ext = BRA.e, main=paste("meanDEF",season_label,"rcp",RCP,"_",year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    #dev.off()
    
    #pdf(paste0(outputDir,"/",className,"/meanPET_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(avPET[["index_1"]], ext = BRA.e, main=paste("meanPET",season_label,"rcp",RCP,"_",year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    #dev.off()
    
    #pdf(paste0(outputDir,"/",className,"/meanDI_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(avDi[["index_1"]], ext = BRA.e, main=paste("meanDI",season_label,"rcp",RCP,"_",year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    #dev.off()
    
    #pdf(paste0(outputDir,"/",className,"/DEFmonths_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(DEFmonths[["index_1"]], ext = BRA.e, main=paste("count DEFmonths",season_label,"rcp",RCP,"_",year,hemi, sep=" "))
    #dev.off()
    
    #pdf(paste0(outputDir,"/",className,"/DEFmonths_prop_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(DEFmonths_prop[["index_1"]], ext = BRA.e, main=paste("prop DEFmonths",season_label,"rcp",RCP,"_",year,hemi, sep=" "))
    dev.off()
    
  }
  
  #create the Moisture Capital Map
  MoistureCap <- round(((75 - avDi[["index_1"]]) / 75),3)
  MoistureCap[MoistureCap[]<0]<-0
  
  
  if(GS) {
    GSCap <- MoistureCap
    GSCap[munis.r %/% 100000 == 42]<-0    #set SC state to 0
    GSCap[munis.r %/% 100000 == 43]<-0    #set RS state to 0
    
    pdf(paste0(outputDir,"/",className,"/GSCap_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(GSCap, ext = BRA.e, main=paste("GSCap",season_label,"rcp",RCP,"_",year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    dev.off()
    
    writeRaster(GSCap, paste0(outputDir,"/",className,"/GSCap_",season_label,"_",hemi,"_rcp",RCP,"_",year,".asc"), format = 'ascii', overwrite=T)
    
  }
  
  if(!GS){
    
    writeRaster(MoistureCap, paste0(outputDir,"/",className,"/MoistureCap_",season_label,"_",hemi,"_rcp",RCP,"_",year,".asc"), format = 'ascii', overwrite=T)
    
    pdf(paste0(outputDir,"/",className,"/MoistureCap_",season_label,"_","rcp",RCP,"_",year,hemi,".pdf"))
    plot(MoistureCap, ext = BRA.e, main=paste("MoistureCap",season_label,"rcp",RCP,"_",year,hemi, sep=" "))  #need to use "index_1" to get to months labelled 1 in season_indices
    dev.off()
  }
  
  rm(pre,tmn,tmx,pre.b,tmn.b,tmx.b,PET.b,DEF.b,ET.b)
  
}


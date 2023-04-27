#summarise gridded ClimBRA daily files to monthly
#for gridded data from https://doi.org/10.57760/sciencedb.02316

library(terra)
library(zoo)  #required for as.yearmon (masks terra::time)
library(lubridate)  #for ym()
library(logger)
library(ncdf4)

path <- "/home/james/Documents/data/"
logfile <- "ClimBRA-summary.log"
log_appender(appender_tee(logfile))

#set model, scenarios, start, end to summarise
model <- c("EC-EARTH3")
scenario <- c("ssp585")
startym <- "2018-01"
endym <- "2018-03"
log_info(model)
log_info(scenario)

for(m in model){
  for(s in scenario){
    
    #precipitation with sum
    fpin <- paste0(path,model,"-pr-",scenario,".nc")
    
    log_info(paste0("reading ",fpin))
    dat <- rast(fpin)
    
    dat <- subset(dat,time(dat) > ym(startym) & time(dat) <= ym(endym))
    
    log_info("starting summary (year-months sum)")
    datym <- tapp(dat, "yearmonths", sum)
    #terra::time(datym) <- seq(as.yearmon(startym),as.yearmon(endym),by=1/12)
    #terra::time(datym) <- seq(ym(startym),ym(endym),by="months")
    #terra::time(datym, tstep="yearmonths") <- seq(ym(startym),ym(endym),by="months")
    #names(datym) <- seq(ym(startym),ym(endym),by="months")
    #names(datym) <- seq(as.yearmon(startym),as.yearmon(endym),by=1/12)
    
    fpout <- paste0(path,model,"-pr-",scenario,"-months-",startym,"-to-",endym,".nc")
    log_info(paste0("writing ", fpout))
    writeCDF(datym, fpout, overwrite=T,
             varname="pr", unit="mm", 
             longname=paste0("total monthly precipitation ",startym," to ",endym))
    
    #terra cannot write year-month to nc as it uses GDAL engine
    #terra also does not preserve layer names (unless they are defined as sub-datasets)
    #also zvals must be values (cannot be times)
    #so best we can do is save as much info in globals (and filename) and then 
  }
}



    #temperature variables with mean
    for(tas in c("tasmax", "tasmin")){
      
      fpin <- paste0(path,model,"-",tas,"-",scenario,".nc")
      log_info(paste0("reading ",fpin))
      dat <- rast(fpin)
      
      dat <- subset(dat,time(dat) > ym(startym) & time(dat) <= ym(endym))
      
      log_info("starting summary (year-months mean)")
      datym <- tapp(dat, "yearmonths", mean)
      terra::time(datym, tstep="yearmon") <- seq(ym(startym),ym(endym),by="months")

      nctime <- ncdim_def("Time","yearmon", seq(ym(startym),ym(endym),by="months"))
      
      fpout <- paste0(path,model,"-",tas,"-",scenario,"-months.nc" )
      log_info(paste0("writing ", fpout))
      writeCDF(datym, fpout, overwrite=T)
      
      zvals = lubridate::parse_date_time(seq(ym(startym),ym(endym),by="months"),
                                         orders = 'ymd', tz = "UTC")
      zvals = as.integer(zvals)
      
      writeCDF(datym, fpout, overwrite=T, 
               varname=tas, unit="degC", 
               longname=paste0("mean monthly ", tas), 
               zname="Time", ncdf4::ncvar_put(nc, 'Time', zvals))
               
    }
  }
}


setwd('/datasets/work/oa-alantis/work/SS_Hydro/connie_atlantis/Hydrodynamic/')
library('ncdf4')
file   <- 'https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DTracerFields1hV19-05.nc'
u.file <- 'https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DuGridFields1hV19-05.nc?uVelocity'
v.file <- 'https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DvGridFields1hV19-05.nc?vVelocity'
w.file <- 'https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSg3DwGridFields1hV19-05.nc?wVelocity'
dimensions <- '[(0.5000003):1:(441.4661)][(0.0):1:(897.0)][(0.0):1:(397.0)]'
date.bgn <- as.Date('2015-01-01')
date.fin <- as.Date('2018-01-01') #as.Date('2019-05-21')
while(date.fin >= (date.bgn + 3)){
  cat(paste('\n Reading day ', date.bgn, '\n'))
  date.end  <- date.bgn + 3
  t.record  <- paste0('[(', date.bgn, 'T00:30:00Z):6:(', date.end, 'T23:30:00Z)]')
  ## Salinity and Temperature
  download.file(url= paste0(file,'?salinity', t.record, dimensions,',temperature', t.record, dimensions), destfile=paste0('Raw_Variables_data/', date.bgn, '_Raw_variables.nc'))
  ## U-velocity
  download.file(url= paste0(u.file, t.record, dimensions), destfile=paste0('Raw_transport_data/', date.bgn, '_URaw_variables.nc'))
  ## V-velocity
  download.file(url= paste0(v.file, t.record, dimensions), destfile=paste0('Raw_transport_data/', date.bgn, '_VRaw_variables.nc'))
  ## W - Variable
  download.file(url= paste0(w.file, t.record, '[(0.0):1:(428.0)][(0.0):1:(897.0)][(0.0):1:(397.0)]'), destfile=paste0('Raw_Variables_data/', date.bgn, '_WRaw_variables.nc'))
  date.bgn <- date.end+1
}


setwd('/datasets/work/oa-alantis/work/SS_Hydro/connie_atlantis/Hydrodynamic/')
library('ncdf4')
file   <- 'https://salishsea.eos.ubc.ca/erddap/griddap/ubcSSaSurfaceAtmosphereFieldsV1.nc'
dimensions <- '[(0.0):1:(662500.0)][(0.0):1:(637500.0)]'
date.bgn <- as.Date('2015-01-01')
date.fin <- as.Date('2018-01-01')
while(date.fin >= (date.bgn + 3)){
  cat(paste('\n Reading day ', date.bgn, '\n'))
  date.end  <- date.bgn + 3
  t.record  <- paste0('[(', date.bgn, 'T00:30:00Z):6:(', date.end, 'T23:30:00Z)]')
  ## Wind
  download.file(url= paste0(file,'?u_wind', t.record, dimensions,',v_wind', t.record, dimensions), destfile=paste0('wind/', date.bgn, '_Wind_variables.nc'))
  date.bgn <- date.end+1
}

% Water layer thickness
addpath('/home/por07g/Documents/PhD/Atlantis_Model/tools/physics/')
layerdepth = [0  25 50 100 250 400 700]; %% This structure is related with
dlev       = [0 diff(layerdepth)]
sum(dlev)
numLayers = get_numLayers('SalishSea_July172019_utm_fix.bgm', dlev)

numLayers
nc=netcdf('JFRE_2000temp_Bec.nc')

ncdump(nc)
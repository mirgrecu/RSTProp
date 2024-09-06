import lkTables
import numpy as np
import glob
lookupT=lkTables.scattTables()

bwy=[32.1,32.1,18.1,18.1,16.,15.6,15.6,7.2,7.2]
bwx=[19.4,19.4,10.9,10.9,9.7,9.4,9.4,4.4,4.4]
bwx=np.array(bwx)
bwy=np.array(bwy)

def antenna_pattern(bwx, bwy):
    bpatt=np.zeros((9,7),float)
    for i in range(9):
        for j in range(7):
            ddx=(i-4)*5
            ddy=(j-3)*5
            y2=(((ddx/bwx)**2+(ddy/bwy)**2)*4*np.log(2.))
            bpatt[i,j]=np.exp(-y2)
    return bpatt

ant_pattL=[]
for i in range(9):
    ant_patt=antenna_pattern(bwx[i], bwy[i])
    ant_pattL.append(ant_patt/ant_patt.sum())

import netCDF4 as nc
def readCMB(fname): # reads relevant data from the CMB file
    fh_cmb=nc.Dataset(fname)
    qv=fh_cmb["KuKaGMI/vaporDensity"][:,:,:]
    press=fh_cmb["KuKaGMI/airPressure"][:,:,:]
    envNodes=fh_cmb["KuKaGMI/envParamNode"][:,:,:]
    airTemp=fh_cmb["KuKaGMI/airTemperature"][:,:,:]
    skTemp=fh_cmb["KuKaGMI/skinTemperature"][:,:]
    binNodes=fh_cmb["KuKaGMI/phaseBinNodes"][:,:]
    pwc=fh_cmb["KuKaGMI/precipTotWaterCont"][:,:,:]
    sfcEmiss=fh_cmb["KuKaGMI/surfEmissivity"][:,:,:]
    dm=fh_cmb["KuKaGMI/precipTotDm"][:,:,:]
    cldw=fh_cmb["KuKaGMI/cloudLiqWaterCont"][:,:,:]
    sfcBin=fh_cmb["KuKaGMI/Input/surfaceRangeBin"][:,:,:]
    zCorrected=fh_cmb["KuGMI/correctedReflectFactor"][:,:,:]
    pType=fh_cmb["KuKaGMI/Input/precipitationType"][:,:]
    lon=fh_cmb["KuKaGMI/Longitude"][:,:]
    lat=fh_cmb["KuKaGMI/Latitude"][:,:]
    return qv,press,envNodes,airTemp,skTemp,binNodes,pwc,sfcEmiss,dm,cldw,sfcBin,zCorrected,pType,lon,lat

def readGMI(fname):
    nc_gmi=nc.Dataset(fname)
    gmi_lat=nc_gmi["S1/Latitude"][:,:]
    gmi_lon=nc_gmi["S1/Longitude"][:,:]
    gmi_tc=nc_gmi["S1/Tc"][:,:]
    return gmi_lat, gmi_lon, gmi_tc

from pyresample import geometry, image, kd_tree
from pyresample.kd_tree import resample_custom

wf = lambda r: 1 - r/20000.0

def resample_gmi(gmi_lat, gmi_lon, gmi_tc, lon, lat):
    # define the GMI grid
    input_def = geometry.SwathDefinition(lons=gmi_lon[:,:], lats=gmi_lat[:,:])
    output_def = geometry.SwathDefinition(lons=lon, lats=lat)
# Resample the tb_s1 data to the CMB grid using gaussian resampling

    gmi_tc_resampled = resample_custom(input_def, gmi_tc[:,:,:], output_def, radius_of_influence=30000, neighbours=10, weight_funcs=[wf for k in range(9)], fill_value=None)
    return gmi_tc_resampled
import radtran as rt
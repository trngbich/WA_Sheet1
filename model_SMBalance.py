# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 08:22:17 2019

@author: sse
"""
import os
import numpy as np
import gdal
import xarray as xr
#import glob
#import datetime
import warnings
import time


#%%
import os
import numpy as np
import gdal
import xarray as xr
import time

#%%
def open_nc(nc,timechunk=1,chunksize=1000):
    dts=xr.open_dataset(nc)
    key=list(dts.keys())[0]
    var=dts[key].chunk({"time": timechunk, "latitude": chunksize, "longitude": chunksize}) #.ffill("time")
    return var,key

def SCS_calc_SRO(P,I,NRD,SMmax,SM, cf):    
    SRO = ((((P-I)/NRD)**2)/((P-I)/NRD+cf*(SMmax-SM))).where(P-I>0,P*0)
    return SRO*NRD

def get_rootdepth(version = '1.0'):
    
    lcc_code = dict()
    lcc_code['1.0'] = {
       'Shrubland':20,
       'Grassland':30,
       'Rainfed Cropland':41,
       'Irrigated Cropland':42,
       'Fallow Cropland':43,
       'Built-up':50,
       'Bare/sparse vegetation':60,
       'Permanent snow/ ice':70,
       'Water bodies':80,
       'Temporary water bodies':81,
       'Shrub or herbaceous cover, flooded':90,     
       'Tree cover: closed, evergreen needle-leaved':111,
       'Tree cover: closed, evergreen broad-leaved':112, 
       'Tree cover: closed, deciduous broad-leaved':114, #
       'Tree cover: closed, mixed type':115, #
       'Tree cover: closed, unknown type':116, #
       'Tree cover: open, evergreen needle-leaved':121,#
       'Tree cover: open, evergreen broad-leaved':122, #
       'Tree cover: open, deciduous needle-leaved':123, #
       'Tree cover: open, deciduous broad-leaved':124, #
       'Tree cover: open, mixed type':125, #
       'Tree cover: open, unknown type':126, #     
       'Seawater':200, #       
       }
    
    root_depth = dict()
    '''
    based on Global estimation of effective plant rooting depth: 
    Implications for hydrological modeling by Yang et al (2016)
    '''
    root_depth['1.0'] = {
       'Shrubland':370,
       'Grassland':510,
       'Rainfed Cropland':550,
       'Irrigated Cropland':550,
       'Fallow Cropland':550,
       'Built-up':370,
       'Bare/sparse vegetation':370,
       'Permanent snow/ ice':0,
       'Water bodies':0,
       'Temporary water bodies':0,
       'Shrub or herbaceous cover, flooded':0,     
       'Tree cover: closed, evergreen needle-leaved':1800,
       'Tree cover: closed, evergreen broad-leaved':3140, 
       'Tree cover: closed, deciduous broad-leaved':1070, #
       'Tree cover: closed, mixed type':2000, #
       'Tree cover: closed, unknown type':2000, #
       'Tree cover: open, evergreen needle-leaved':1800,#
       'Tree cover: open, evergreen broad-leaved':3140, #
       'Tree cover: open, deciduous needle-leaved':1070, #
       'Tree cover: open, deciduous broad-leaved':1070, #
       'Tree cover: open, mixed type':2000, #
       'Tree cover: open, unknown type':2000, #     
       'Seawater':0, #       
    }
    
    return lcc_code[version], root_depth[version]


def root_dpeth(lu):
    rootdepth=lu.copy()
    rootdepth.name='Root depth'
    rootdepth.attrs={'units':'mm',
                    'quantity':'Effective root depth',
                    'source':'Root depth lookup table',
                    'period':'year'}
    lu_categories, root_depth = get_rootdepth(version = '1.0')
    for key in root_depth.keys():
      lu_code=lu_categories[key]
      rd = root_depth[key]
      rootdepth=rootdepth.where(lu!=lu_code,rd)
    return rootdepth 


def get_fractions(version = '1.0'):
    consumed_fractions = dict()
    
    consumed_fractions['1.0'] = {
       'Shrubland':1.00,
       'Grassland':1.00,
       'Rainfed Cropland':1.00,
       'Irrigated Cropland':0.80,
       'Fallow Cropland':1.00,
       'Built-up':1.00,
       'Bare/sparse vegetation':1.00,
       'Permanent snow/ ice':1.00,
       'Water bodies':1.00,
       'Temporary water bodies':1.00,
       'Shrub or herbaceous cover, flooded':1.00,     
       'Tree cover: closed, evergreen needle-leaved':1.00,
       'Tree cover: closed, evergreen broad-leaved':1.00, 
       'Tree cover: closed, deciduous broad-leaved':1.00, #
       'Tree cover: closed, mixed type':1.00, #
       'Tree cover: closed, unknown type':1.00, #
       'Tree cover: open, evergreen needle-leaved':1.00,#
       'Tree cover: open, evergreen broad-leaved':1.00, #
       'Tree cover: open, deciduous needle-leaved':1.00, #
       'Tree cover: open, deciduous broad-leaved':1.00, #
       'Tree cover: open, mixed type':1.00, #
       'Tree cover: open, unknown type':1.00, #     
       'Seawater':1.00, #       

    }
    
    return consumed_fractions[version]

def Consumed_fraction(lu):
    f_consumed=lu.copy()
    f_consumed.name='Consumed fraction'
    f_consumed.attrs={'units':'Fraction',
                    'quantity':'Consumed fraction',
                    'source':'Consumed fraction look-up table',
                    'period':'year'}
    consumed_fractions = get_fractions(version = '1.0')
    lu_categories, root_depth = get_rootdepth(version = '1.0')
    for key in consumed_fractions.keys():
        lu_code=lu_categories[key]
        consumed_fraction = consumed_fractions[key]
        f_consumed = f_consumed.where(lu!=lu_code,consumed_fraction)
    return f_consumed 

def OpenAsArray(fh, bandnumber = 1, dtype = 'float32', nan_values = False):
    """
    Open a map as an numpy array. 
    
    Parameters
    ----------
    fh: str
        Filehandle to map to open.
    bandnumber : int, optional 
        Band or layer to open as array, default is 1.
    dtype : str, optional
        Datatype of output array, default is 'float32'.
    nan_values : boolean, optional
        Convert he no-data-values into np.nan values, note that dtype needs to
        be a float if True. Default is False.
        
    Returns
    -------
    Array : ndarray
        Array with the pixel values.
    """
    datatypes = {"uint8": np.uint8, "int8": np.int8, "uint16": np.uint16, "int16":  np.int16, "Int16":  np.int16, "uint32": np.uint32,
    "int32": np.int32, "float32": np.float32, "float64": np.float64, "complex64": np.complex64, "complex128": np.complex128,
    "Int32": np.int32, "Float32": np.float32, "Float64": np.float64, "Complex64": np.complex64, "Complex128": np.complex128,}
    DataSet = gdal.Open(fh, gdal.GA_ReadOnly)
    Type = DataSet.GetDriver().ShortName
    if Type == 'HDF4':
        Subdataset = gdal.Open(DataSet.GetSubDatasets()[bandnumber][0])
        NDV = int(Subdataset.GetMetadata()['_FillValue'])
    else:
        Subdataset = DataSet.GetRasterBand(bandnumber)
        NDV = Subdataset.GetNoDataValue()
    Array = Subdataset.ReadAsArray().astype(datatypes[dtype])
    if nan_values:
        Array[Array == NDV] = np.nan
    Array = Array.astype(np.float32)
    return Array


#%% main
def run_SMBalance(MAIN_FOLDER,p_in,e_in,i_in,rd_in,lu_in,smsat_file,
         start_year=2009,f_perc=1,f_Smax=0.9, cf =  20,
         chunks=[1,1000,1000]):
    '''
    Arguments:
        
    ## required   
    MAIN_FOLDER='$PATH/nc/'
    p_in = '$PATH/p_monthly.nc'
    e_in = '$PATH/e_monthly.nc'
    i_in = '$PATH/i_monthly.nc'
    rd_in = '$PATH/nRD_monthly.nc'
    lu_in = '$PATH/lcc_yearly.nc'
    smsat_file = '$PATH/p_monthly.nc'
    start_year=2009
    
    #default
    f_perc=1 # percolation factor
    f_Smax=0.9 #threshold for percolation
    cf =  20 #f_Ssat soil mositure correction factor to componsate the variation in filling up and drying in a month
 
    '''
    warnings.filterwarnings("ignore", message='invalid value encountered in greater')
    warnings.filterwarnings("ignore", message='divide by zero encountered in true_divide')
    warnings.filterwarnings("ignore", message='invalid value encountered in true_divide')
    warnings.filterwarnings("ignore", message='overflow encountered in exp')

    tchunk=chunks[0]
    chunk=chunks[1]
    Pt,key=open_nc(p_in,timechunk=tchunk,chunksize=chunk)
    E,key=open_nc(e_in,timechunk=tchunk,chunksize=chunk)
    Int,key=open_nc(i_in,timechunk=tchunk,chunksize=chunk)
    nRD,key=open_nc(rd_in,timechunk=tchunk,chunksize=chunk)
    LU,key=open_nc(lu_in,timechunk=tchunk,chunksize=chunk)
#    thetasat=xr.open_dataset(thetasat_in).thetasat.chunk({'latitude':chunk,'longitude':chunk})
    ### Get thetasat DataArray from tiff
    thetasat_data=OpenAsArray(smsat_file,nan_values=True)
    thetasat=xr.DataArray(thetasat_data, coords = {"latitude": Pt.latitude,
                                                    "longitude":Pt.longitude}, 
                           dims = ["latitude", "longitude"])
    thetasat=thetasat.chunk({"latitude": chunk, "longitude": chunk}) #.ffill("time")
    
    # del thetasat_data
    
    ### convert nRD = 0 to 1
    nRD = nRD.where(nRD!=0,1)
    
    ### Create 
    SM=E[0]*0
    
    for j in range(len(LU.time)): 
        t1 = j*12
        t2 = (j+1)*12    
       
        lu = LU.isel(time=j)
        
        #mask lu for water bodies
        mask = xr.where(((lu==80) | (lu==81) | (lu==70) | (lu==200)), 1,np.nan)
        
        Rd = root_dpeth(lu)
        SMmax=thetasat*Rd    
        f_consumed = Consumed_fraction(lu)    
        for t in range(t1,t2):
            print('time: ', t)
            SMt_1=SM 
            P = Pt.isel(time=t)
            ETa = E.isel(time=t)
            I = Int.isel(time=t)
            NRD = nRD.isel(time=t)
            
            ### calculate surface runoff  as a function of SMt_1
            SRO=SCS_calc_SRO(P,I,NRD,SMmax,SMt_1,cf)        
             
            ### calculate Percolation as a function of SMt_1
            perc=(SMt_1*(xr.ufuncs.exp(-f_perc/SMt_1))).where(SMt_1>f_Smax*SMmax,P*0)
                    
            ### Calculate SM temp
            Stemp = SMt_1+(P-I)-(ETa-I)-SRO-perc
            
            
            ### Calculate ETincr, ETrain, Qsupply, and update SM
            ETincr = (P*0).where(Stemp>=0, -1*Stemp)
           
            # adjust ETincr for water bodies
            ETincr = (P*0).where(((mask==1) & (P-ETa >= 0)), 
                      (ETa-P).where(((mask==1) & (P-ETa <0)), ETincr))
            
            ETrain = ETa.where(Stemp>=0,ETa-ETincr)
            Qsupply = (P*0).where(Stemp>=0,ETincr/f_consumed)
            SM = Stemp.where(Stemp>=0,Stemp+Qsupply)
            
            ### Calculate increametal percolation and increamental runoff
            perc_incr = (SM-SMmax)*perc/(perc+SRO).where(SM>SMmax, P*0)
            SROincr = (SM-SMmax-perc_incr).where(SM>SMmax, P*0)
    #        overflow = SM-SMmax # we don't use overflow at 
            SM=SM.where(SM<SMmax, SMmax)
    #        SRO=SRO+overflow.where(overflow>0, SRO)
          
    
    
            if t == 0:
                etb = ETincr
                etg = ETrain
            else:
                etb = xr.concat([etb, ETincr], dim='time')  
                etg = xr.concat([etg, ETrain], dim='time')
            attrs={"units":"mm/month", "source": "FAO WaPOR", "quantity":"Rainfall_ET_M"}
            etg.attrs=attrs
            etg.name = 'Rainfall_ET_M'
    
            attrs={"units":"mm/month", "source": "FAO WaPOR", "quantity":"Increamental_ET_M"}
            etb.attrs=attrs
            etb.name = 'Increamental_ET_M'
    
            ###
#            timechunk = 1
#            chunks = [timechunk,latchunk, lonchunk]
            comp = dict(zlib=True, complevel=9, least_significant_digit=2, chunksizes=chunks)
    
    
            ##green ET 
            print("\n\nwriting the ET_rain netcdf file\n\n")
            etg_outfolder=os.path.join(MAIN_FOLDER,'etg')
            if not os.path.exists(etg_outfolder):
              os.makedirs(etg_outfolder)
            nc_fn=r'et_g_monthly_4.nc'
            nc_path=os.path.join(etg_outfolder,nc_fn)
            encoding = {"Rainfall_ET_M": comp}
            etg.to_netcdf(nc_path,encoding=encoding)    

        del etg
        #
        ###Blue ET datase
        print("\n\nwriting the ET_increamental netcdf file\n\n")
        etb_outfolder=os.path.join(MAIN_FOLDER,'etb')
        if not os.path.exists(etb_outfolder):
          os.makedirs(etb_outfolder)
        nc_fn=r'et_b_monthly_4.nc'
        nc_path=os.path.join(etb_outfolder,nc_fn)
        encoding = {"Increamental_ET_M": comp}
        etb.to_netcdf(nc_path,encoding=encoding)

        del etb
        
        

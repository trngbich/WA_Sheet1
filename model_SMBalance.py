# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 08:22:17 2019

@author: sse
"""
import os
import numpy as np
import gdal
import xarray as xr
import glob
import datetime
import warnings


#%% functions
def open_nc(nc,timechunk=1,chunksize=1000):
    dts=xr.open_dataset(nc)
    key=list(dts.keys())[0]
    var=dts[key].chunk({"time": timechunk, "latitude": chunksize, "longitude": chunksize}) #.ffill("time")
    return var,key

def calc_SM_bucket(P,ETa,I,SMmax,SMgt_1,SMincrt_1,f_consumed):
    SMt_1=SMgt_1+SMincrt_1
    ETr=np.where(I>0,ETa-I, ETa)
    SMtemp=SMgt_1+np.maximum(P-I,P*0)
    SMg=np.where(SMtemp-ETr>0,SMtemp-ETr,P*0)
    ETincr=np.where(SMtemp-ETr>0,P*0,ETr-SMtemp)
    Qsupply=np.where(ETincr>SMincrt_1,(ETincr-SMincrt_1)/f_consumed, P*0) 
            
    #SMincr=np.where(ETincr>SMincrt_1,SMincrt_1+Qsupply-ETincr,SMincrt_1-ETincr)
    SMincr=SMincrt_1+Qsupply-ETincr
    SMincr=np.where(SMg+SMincr>SMmax,SMincr-SMincr/(SMg+SMincr)*((SMg+SMincr)-SMmax),SMincr)     
    SMg=np.where(SMg+SMincr>SMmax,(SMg-SMg/(SMg+SMincr)*((SMg+SMincr)-SMmax)),SMg)
    
    #SMincr=np.where(SMg+SMincr>SMmax,SMincr-((SMg+SMincr)-SMmax),SMincr)
    SM=np.where(SMg+SMincr>SMmax,SMmax,SMg+SMincr)
    dsm=SM-SMt_1
    ETg=ETa-ETincr
    return SMg,SMincr,SM,dsm,Qsupply,ETincr,ETg

def calc_SCS_SRO(P,I,SMmax,SM,Qsupply, cf):    
    '''
    SCS formula
    '''
    SM=np.where(SM>SMmax,SMmax,SM)        
    SRO= np.where((P-I+Qsupply)>0,((P-I+Qsupply)**2)/(P-I+Qsupply+cf*(SMmax-SM)),P*0)
    SROg= np.where(P-I>0,((P-I)**2)/(P-I+cf*(SMmax-SM)),P*0)
    SROincr=SRO-SROg
    return SRO,SROincr

def get_rootdepth(version = '1.0'):
    '''
    rootdepth lookup table
    '''
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

def calc_rootdepth(lu):
    rootdepth=np.copy(lu)
    lu_categories, root_depth = get_rootdepth(version = '1.0')
    for key in root_depth.keys():
        mask = np.logical_or.reduce([lu == lu_categories[key]])
        rd = root_depth[key]
        rootdepth[mask] = rd
    return rootdepth 


def get_fractions(version = '1.0'):
    '''
    consumed fraction lookup table
    '''
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
#            'Forests':              1.00,
#            'Shrubland':            1.00,
#            'Rainfed Crops':        1.00,
#            'Forest Plantations':   1.00,
#            'Natural Water Bodies': 1.00,
#            'Wetlands':             1.00,
#            'Natural Grasslands':   1.00,
#            'Other (Non-Manmade)':  1.00,#0.4
#            'Irrigated crops':      0.80,#0.8
#            'Managed water bodies': 1.00,#0.4
#            'Other':                1.00,#0.4
#            'Residential':          1.00,
#            'Greenhouses':          0.95,
#            'Aquaculture':          1.00
    }
    
    return consumed_fractions[version]

def calc_consumed_fraction(lu):
    f_consumed=np.copy(lu)
    consumed_fractions = get_fractions(version = '1.0')
    lu_categories, root_depth = get_rootdepth(version = '1.0')
    for key in consumed_fractions.keys():
        mask = np.logical_or.reduce([lu == lu_categories[key]])
        consumed_fraction = consumed_fractions[key]
        f_consumed[mask] = consumed_fraction
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
    return Array


#%% main
def run_SMBalance(MAIN_FOLDER,p_in,e_in,i_in,rd_in,lu_in,smsat_file,
         start_year=2009,f_perc=1,f_Smax=0.9, cf =  20):
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
    
    start_time=datetime.datetime.now()
    
    Pt,_=open_nc(p_in,chunksize=1000)
    E,_=open_nc(e_in,chunksize=1000)
    Int,_=open_nc(i_in,chunksize=1000)
    nRD,_=open_nc(rd_in,chunksize=1000)
    LU,_=open_nc(lu_in,chunksize=1000)
    
    ##read saturation file
    thetasat = OpenAsArray(smsat_file,nan_values=True)    
    nrd = nRD.where(nRD!=0,1)
     
    SMg=E[0]*0
    SMincr=E[0]*0
    
    for j in range(len(LU.time)): 
        t1 = j*12
        t2 = (j+1)*12    
        etb = np.zeros((12,E[0].shape[0],E[0].shape[1]))
        etg = np.zeros((12,E[0].shape[0],E[0].shape[1]))
        lu = LU.isel(time=j).values
        Rd = calc_rootdepth(lu)
        SMmax=thetasat*Rd    
        f_consumed = calc_consumed_fraction(lu)    
        for t in range(t1,t2):
            print('time: ', t)
            if t==0:
                SMgt_1=E[0]*0
                Qsupply=E[0]*0
                SMincrt_1=E[0]*0
            else:
                SMgt_1=SMg 
                SMincrt_1=SMincr 
                
            p = Pt.isel(time=t).values
            e = E.isel(time=t).values
            i = Int.isel(time=t).values
            NRD = nrd.isel(time=t).values
            P = p/NRD
            ETa = e/NRD
            I = i/NRD
            
            ### Step 1: Soil moisture
            SMg,SMincr,SM,dsm,Qsupply,ETincr,ETg=calc_SM_bucket(P,ETa,I,SMmax,
                                                                SMgt_1,SMincrt_1,
                                                                f_consumed) 
            
            SRO,SROincr=calc_SCS_SRO(P,I,SMmax,SM,Qsupply,cf)        
            ### Step 3: Percolation       
           
                   
            perc=np.where(SM>f_Smax*SMmax,SM*(np.exp(-f_perc/SM)),P*0)
            SMincr = np.where(SMincr-SROincr*NRD>0,SMincr-SROincr*NRD,P*0)
            SMg_ratio = np.where(SM==0,1,SMg/SM)
            SMincr_ratio = np.where(SM==0,1,SMincr/SM)
            perc_green = perc*SMg_ratio
            perc_incr = perc*SMincr_ratio
            perc_green = np.where(SMg>perc_green, perc_green, SMg*0)
            perc_incr = np.where(SMincr>perc*SMincr_ratio, perc_incr, SMincr*0)
            
            SMg = np.where(SMg-perc_green>0, SMg-perc_green, P*0)
            SMincr = np.where(SMincr-perc_incr >0, SMincr-perc_incr, P*0)
            
            ###updating the soil moisture by subtracting the surface runoff
            SROg=SRO-SROincr
            SMg=np.where(SMg-SROg>0,SMg-SROg,SMg*0)
            SMincr=np.where(SMincr-SROincr>0,SMincr-SROincr,SMincr*0)
    
            SM = SMg+SMincr
            
            ### Step 5: Store monthly data of the year
            k = int(t-(j*12))
            etb[k,:,:] = ETincr*NRD
            etg[k,:,:] = ETg*NRD
            
            del p
            del ETa
            del I
            del NRD
            del SMgt_1 
            del Qsupply 
            del SMincrt_1
            
        year = start_year+j    
        time_ds = E.time[t1:t2]
        ds = xr.Dataset({})
        ds['etincr'] = (('time','latitude', 'longitude'), etb)
        ds['etrain'] = (('time','latitude', 'longitude'), etg)
        ds = ds.assign_coords(time = time_ds, latitude=E.latitude,longitude=E.longitude )
        
        ##rainfall ET data
        et_rain = ds.etrain
        attrs={"units":"mm/month", "source": "SM Balance model", 
               "quantity":"ET rainfall"}
        et_rain.assign_attrs(attrs)
        et_rain.name = "etrain"
        et_rain_dts=et_rain.chunk({"latitude":-1,"longitude":-1}).to_dataset()
        etrain_outfolder=os.path.join(MAIN_FOLDER,'etrain')
        if not os.path.exists(etrain_outfolder):
          os.makedirs(etrain_outfolder)
        nc_fn=r'et_rain_monthly_'+str(year)+'.nc'
        nc_path=os.path.join(etrain_outfolder,nc_fn)
        comp = dict(zlib=True, complevel=9, least_significant_digit=3)
        encoding = {"etrain": comp}
        et_rain_dts.to_netcdf(nc_path,encoding=encoding)
        
        ##Incremental ET data
        et_incr = ds.etincr
        attrs={"units":"mm/month", "source": "SM Balance model", 
               "quantity":"ET incremental"}
        et_incr.assign_attrs(attrs)
        et_incr.name = "etincr"
        et_incr_dts=et_incr.chunk({"latitude":-1,"longitude":-1}).to_dataset()
        etincr_outfolder=os.path.join(MAIN_FOLDER,'etincr')
        if not os.path.exists(etincr_outfolder):
          os.makedirs(etincr_outfolder)
        nc_fn=r'et_incr_monthly_'+str(year)+'.nc'
        nc_path=os.path.join(etincr_outfolder,nc_fn)
        et_incr_dts.to_netcdf(nc_path,encoding={"etincr":{'zlib':True}})
        
        del etg
        del etb

    elapsed_time=datetime.datetime.now()-start_time
    print('Model finished in ',elapsed_time)        
    print('Output1 ETrain folder: ',etrain_outfolder)
    print('Output2 ETincr folder: ',etincr_outfolder)
    return etrain_outfolder,etincr_outfolder
        
def merge_yearly_nc(nc_folder,out_nc,varname=None):    
    fhs=glob.glob(os.path.join(nc_folder,'*.nc'))
    arr_list=[]
    for fh in fhs:
        arr,key=open_nc(fh)
        arr_list.append(arr)
    var=xr.concat(arr_list, dim='time')
    if varname is None:
        varname=key
    attrs={"units":"mm/month", "source": "Merged yearly netCDF files    from {0}".format(nc_folder), "quantity":varname}
    var.assign_attrs(attrs)
    dts=var.chunk({"latitude":-1,"longitude":-1}).to_dataset()  
    comp = dict(zlib=True, complevel=9, least_significant_digit=3)
    encoding = {key: comp}
    dts.to_netcdf(out_nc,encoding=encoding)
    print('Finish merging ',out_nc)
    
        

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
def open_nc(nc,timechunk=-1,chunksize=1000):
    dts=xr.open_dataset(nc)
    key=list(dts.keys())[0]
    var=dts[key].chunk({"time": timechunk, "latitude": chunksize, "longitude": chunksize}) #.ffill("time")
    return var,key

# def SM_bucket(P,ETa,I,SMmax,SMgt_1,SMincrt_1,f_consumed):
#     SMt_1=SMgt_1+SMincrt_1
#     ETr=np.where(I>0,ETa-I, ETa)
#     SMtemp=SMgt_1+np.maximum(P-I,P*0)
#     SMg=np.where(SMtemp-ETr>0,SMtemp-ETr,P*0)
#     ETincr=np.where(SMtemp-ETr>0,P*0,ETr-SMtemp)
#     Qsupply=np.where(ETincr>SMincrt_1,(ETincr-SMincrt_1)/f_consumed, P*0) 
            
#     #SMincr=np.where(ETincr>SMincrt_1,SMincrt_1+Qsupply-ETincr,SMincrt_1-ETincr)
#     SMincr=SMincrt_1+Qsupply-ETincr
#     SMincr=np.where(SMg+SMincr>SMmax,SMincr-SMincr/(SMg+SMincr)*((SMg+SMincr)-SMmax),SMincr)     
#     SMg=np.where(SMg+SMincr>SMmax,(SMg-SMg/(SMg+SMincr)*((SMg+SMincr)-SMmax)),SMg)
    
#     #SMincr=np.where(SMg+SMincr>SMmax,SMincr-((SMg+SMincr)-SMmax),SMincr)
#     SM=np.where(SMg+SMincr>SMmax,SMmax,SMg+SMincr)
#     dsm=SM-SMt_1
#     ETg=ETa-ETincr
#     return SMg,SMincr,SM,dsm,Qsupply,ETincr,ETg
    
def SM_bucket(P,ETa,I,SMmax,SMgt_1,SMincrt_1,f_consumed):
    SMt_1=SMgt_1+SMincrt_1
    ETr=(ETa-I).where(I>0,ETa) 
    SMtemp=SMgt_1+xr.ufuncs.maximum(P-I,P*0) #or SMgt_1+np.maximum(P-I,P*0)
    SMg=(SMtemp-ETr).where(SMtemp-ETr>0,P*0) 
    ETincr= (P*0).where(SMtemp-ETr>0,ETr-SMtemp) 
    Qsupply=((ETincr-SMincrt_1)/f_consumed).where(ETincr>SMincrt_1,P*0)        
    SMincr=SMincrt_1+Qsupply-ETincr
    SMincr=(SMincr-SMincr/(SMg+SMincr)*((SMg+SMincr)-SMmax)).where(SMg+SMincr>SMmax,SMincr)     
    SMg=((SMg-SMg/(SMg+SMincr)*((SMg+SMincr)-SMmax))).where(SMg+SMincr>SMmax,SMg)
    SM=SMmax.where(SMg+SMincr>SMmax,SMg+SMincr)
    dsm=SM-SMt_1
    ETg=ETa-ETincr
    return SMg,SMincr,SM,dsm,Qsupply,ETincr,ETg

# def SCS_calc_SRO(P,I,SMmax,SM,Qsupply, cf):    
#     SM=np.where(SM>SMmax,SMmax,SM)        
#     SRO= np.where((P-I+Qsupply)>0,((P-I+Qsupply)**2)/(P-I+Qsupply+cf*(SMmax-SM)),P*0)
#     SROg= np.where(P-I>0,((P-I)**2)/(P-I+cf*(SMmax-SM)),P*0)
#     SROincr=SRO-SROg
#     return SRO,SROincr

def SCS_calc_SRO(P,I,SMmax,SM,Qsupply, cf):    
    SM=SMmax.where(SM>SMmax,SM)        
    SRO= (((P-I+Qsupply)**2)/(P-I+Qsupply+cf*(SMmax-SM))).where((P-I+Qsupply)>0,P*0)
    SROg= (((P-I)**2)/(P-I+cf*(SMmax-SM))).where(P-I>0,P*0)
    SROincr=SRO-SROg
    return SRO,SROincr

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

# def root_dpeth(lu):
#     rootdepth=np.copy(lu)
#     lu_categories, root_depth = get_rootdepth(version = '1.0')
#     for key in root_depth.keys():
#         mask = np.logical_or.reduce([lu == lu_categories[key]])
#         rd = root_depth[key]
#         rootdepth[mask] = rd
#     return rootdepth 

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

# def Consumed_fraction(lu):
#     f_consumed=np.copy(lu)
#     consumed_fractions = get_fractions(version = '1.0')
#     lu_categories, root_depth = get_rootdepth(version = '1.0')
#     for key in consumed_fractions.keys():
#         mask = np.logical_or.reduce([lu == lu_categories[key]])
#         consumed_fraction = consumed_fractions[key]
#         f_consumed[mask] = consumed_fraction
#     return f_consumed 

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
    return Array


#%% main
def run_SMBalance(MAIN_FOLDER,p_in,e_in,i_in,rd_in,lu_in,smsat_file,
         start_year=2009,f_perc=1,f_Smax=0.9, cf =  20,
         chunk=[-1,1000]):
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
    
    start_time=datetime.datetime.now()
    
    Pt,_=open_nc(p_in)
    E,_=open_nc(e_in,timechunk=chunk[0],chunksize=chunk[1])
    Int,_=open_nc(i_in,timechunk=chunk[0],chunksize=chunk[1])
    nRD,_=open_nc(rd_in,timechunk=chunk[0],chunksize=chunk[1])
    LU,_=open_nc(lu_in,timechunk=chunk[0],chunksize=chunk[1])
    
    ##read saturation file
    thetasat_data=OpenAsArray(smsat_file,nan_values=True)
    thetasat=xr.DataArray(thetasat_data, coords = {"latitude": Pt.latitude,
                                                   "longitude":Pt.longitude}, 
                          dims = ["latitude", "longitude"])
    del thetasat_data
    
    nrd = nRD.where(nRD!=0,1)
    
    SMg=E[0]*0
    SMincr=E[0]*0
    etb = E*0
    etg = E*0
    arr_list_etb=[]
    arr_list_etg=[]
    for j in range(len(LU.time)): 
        t1 = j*12
        t2 = (j+1)*12    
       
        lu = LU.isel(time=j)
        Rd = root_dpeth(lu)
        SMmax=thetasat*Rd    
        f_consumed = Consumed_fraction(lu)    
        for t in range(t1,t2):
            print('time: ', t)
            if t==0:
                SMgt_1=E[0]*0
                Qsupply=E[0]*0
                SMincrt_1=E[0]*0
            else:
                SMgt_1=SMg 
                SMincrt_1=SMincr 
             
             
            p = Pt.isel(time=t)
            e = E.isel(time=t)
            i = Int.isel(time=t)
            NRD = nrd.isel(time=t)
            
            P = p/NRD
            ETa = e/NRD
            I = i/NRD
            
            ### Step 1: Soil moisture
            SMg,SMincr,SM,dsm,Qsupply,ETincr,ETg=SM_bucket(P,ETa,I,SMmax,SMgt_1,SMincrt_1,f_consumed) 
            
            SRO,SROincr=SCS_calc_SRO(P,I,SMmax,SM,Qsupply,cf)        
            ### Step 3: Percolation       
           
                   
            perc=(SM*(np.exp(-f_perc/SM))).where(SM>f_Smax*SMmax,P*0)
            SMincr = (SMincr-SROincr*NRD).where(SMincr-SROincr*NRD>0,P*0)
            SMg_ratio = (SM*0+1).where(SM==0,SMg/SM)
            SMincr_ratio = (SM*0+1).where(SM==0,SMincr/SM)
            perc_green = perc*SMg_ratio
            perc_incr = perc*SMincr_ratio
            perc_green = perc_green.where(SMg>perc_green, SMg*0)
            perc_incr = perc_incr.where(SMincr>perc*SMincr_ratio, SMincr*0)
            
            SMg = (SMg-perc_green).where(SMg-perc_green>0, P*0)
            SMincr = (SMincr-perc_incr).where(SMincr-perc_incr >0,P*0)
            
            ###updating the soil moisture by subtracting the surface runoff
            SROg=SRO-SROincr
            SMg=(SMg-SROg).where(SMg-SROg>0,SMg*0)
            SMincr=(SMincr-SROincr).where(SMincr-SROincr>0,SMincr*0)
    
            SM = SMg+SMincr
            
            ### Step 5: Store monthly data of the year
#            k = int(t-(j*12))
            arr_list_etb.append(ETincr*NRD)
            arr_list_etg.append(ETg*NRD)
            
    etb=xr.concat(arr_list_etb, dim='time')  
    etg=xr.concat(arr_list_etg, dim='time')    
    #year = start_year+j    
    #time_ds = E.time[t1:t2]
    #ds = xr.Dataset({})
    #ds['etb'] = (('time','latitude', 'longitude'), etb)
    #ds['etg'] = (('time','latitude', 'longitude'), etg)
    #ds = ds.assign_coords(time = E.time, latitude=E.latitude,longitude=E.longitude )
    
    
    ###green ET datase
    #et_g = ds.etg
    attrs={"units":"mm/month", "source": "FAO WaPOR", "quantity":"Rainfall_ET_M"}
    etg.assign_attrs(attrs)
    etg.name = "Rainfall_ET_M"
    etg_dts=etg.chunk({"latitude":-1,"longitude":-1}).to_dataset()
    nc_fn=r'et_g_monthly.nc'
    nc_path=os.path.join(MAIN_FOLDER,nc_fn)
    comp = dict(zlib=True, complevel=9, least_significant_digit=3)
    encoding = {"Rainfall_ET_M": comp}
    etg_dts.to_netcdf(nc_path,encoding=encoding)
    #

    ###Blue ET datase
    #et_b = ds.etb
    attrs={"units":"mm/month", "source": "FAO WaPOR", "quantity":"Increamental_ET_M"}
    etb.assign_attrs(attrs)
    etb.name = "Increamental_ET_M"
    etb_dts=etb.chunk({"latitude":-1,"longitude":-1}).to_dataset()
    nc_fn=r'et_b_monthly.nc'
    nc_path=os.path.join(MAIN_FOLDER,nc_fn)
    etb_dts.to_netcdf(nc_path,encoding={"Increamental_ET_M":{'zlib':True}})
    #
    #del etb
    #del etg
        
#def merge_yearly_nc(nc_folder,out_nc,varname=None):     
#    fhs=sorted(glob.glob(os.path.join(nc_folder,'*.nc')))
#    arr_list=[]
#    for fh in fhs:
#        arr,key=open_nc(fh)
#        arr_list.append(arr)
#    var=xr.concat(arr_list, dim='time')
#    if varname is None:
#        varname=key
#    attrs={"units":"mm/month", "source": "Merged yearly netCDF files    from {0}".format(nc_folder), "quantity":varname}
#    var.assign_attrs(attrs)
#    dts=var.chunk({"latitude":-1,"longitude":-1}).to_dataset()  
#    comp = dict(zlib=True, complevel=9, least_significant_digit=3)
#    encoding = {key: comp}
#    dts.to_netcdf(out_nc,encoding=encoding)
#    print('Finish merging ',out_nc)
#    
        

# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 13:42:44 2019

@author: ntr002
WAPOR LCC to LUWA
"""
#os.chdir(r'D:\FAO\WA_Sheet1\WAPORWA') 

import WaporIHE.download.GIS_functions as gis
import numpy as np
import os
import ogr 
import gdal


def Rasterize_shapefile(InputVector, RefImage, OutputImage):
    gdalformat = 'GTiff'
    datatype = gdal.GDT_Byte
    burnVal = 1 #value for the output image pixels
    
    # Get projection info from reference image
    Image = gdal.Open(RefImage, gdal.GA_ReadOnly)
    
    # Open Shapefile
    Shapefile = ogr.Open(InputVector)
    Shapefile_layer = Shapefile.GetLayer()
    
    # Rasterise
    #print("Rasterising shapefile...")
    Output = gdal.GetDriverByName(gdalformat).Create(OutputImage, Image.RasterXSize, Image.RasterYSize, 1, datatype, options=['COMPRESS=DEFLATE'])
    Output.SetProjection(Image.GetProjectionRef())
    Output.SetGeoTransform(Image.GetGeoTransform()) 
    
    # Write data to band 1
    Band = Output.GetRasterBand(1)
    Band.SetNoDataValue(0)
    gdal.RasterizeLayer(Output, [1], Shapefile_layer, burn_values=[burnVal])
    
    # Close datasets
    Band = None
    Output = None
    Image = None
    Shapefile = None

   
def Adjust_GRaND_reservoir(output_raster,WaPOR_LCC,GRaND_Reservoir,
                           Resrv_to_Lake,Lake_to_Reserv):   
  
     #Getting GeoTranformation from LCC map
    driver,NDV,xsize,ysize,GeoT,Projection=gis.GetGeoInfo(WaPOR_LCC)
    #Rasterize selected area for reservoir and un-reservoir shapefile
    Basin_reservoir=os.path.join(os.path.split(Resrv_to_Lake)[0],'Reservoir_GRanD.tif')
    Rasterize_shapefile(GRaND_Reservoir,WaPOR_LCC,Basin_reservoir)    
    Rasterize_shapefile(Resrv_to_Lake,WaPOR_LCC,
                        Resrv_to_Lake.replace('.shp','.tif'))    
    Rasterize_shapefile(Lake_to_Reserv,WaPOR_LCC,
                        Lake_to_Reserv.replace('.shp','.tif'))
    
    #Edit Resvr
    Resrv=gis.OpenAsArray(Basin_reservoir,nan_values=True)
    LCC=gis.OpenAsArray(WaPOR_LCC,nan_values=True)
    UnResrv=gis.OpenAsArray(Resrv_to_Lake.replace('.shp','.tif'),nan_values=True)
    MakeResrv=gis.OpenAsArray(Lake_to_Reserv.replace('.shp','.tif'),nan_values=True)
    
    Resrv=np.where(((LCC==80)*(MakeResrv==1)),1,Resrv)
    Resrv=np.where(((Resrv==1)*(UnResrv==1)),np.nan,Resrv)
    
#    output=os.path.join(os.path.split(Resrv_to_Lake)[0],'Reservoir_adjusted.tif')
    gis.CreateGeoTiff(output_raster, Resrv, driver, NDV, xsize, ysize, GeoT, Projection)
    return output_raster

def Reclass_LCC_to_LUWA(WaPOR_LCC,Output_dir,ProtectedArea_tif,
                        Reservoir_tif,LCC_LUWA_dict=None):    
    if LCC_LUWA_dict is None:
        LCC_LUWA_dict={
                'PLU':(1,[]),
               'ULU':(2,[]),
               'MLU':(3,[41,42]),#Rainfed crop
               'MWU':(4,[42,50]), #irrigated crop and built-up
               } 
    driver,NDV,xsize,ysize,GeoT,Projection=gis.GetGeoInfo(WaPOR_LCC)
    LCC=gis.OpenAsArray(WaPOR_LCC,nan_values=True)
    #ULU: The default is ULU 
    LUWA=2*np.ones(np.shape(LCC),dtype=np.float32)    
    #PLU: WDPA 
    PLU=gis.OpenAsArray(ProtectedArea_tif,nan_values=True)
    LUWA=np.where(PLU==1,1,LUWA)
    #MLU: Rainfed crop => Modified Land Use
    for code in LCC_LUWA_dict['MLU'][1]:
        LUWA=np.where(LCC==code,3,LUWA)
    #MWU: Irrigated crop, Reservoir, Urban => Managed Water Use
    for code in LCC_LUWA_dict['MWU'][1]:
        LUWA=np.where(LCC==code,4,LUWA)
    MWU=gis.OpenAsArray(Reservoir_tif,nan_values=True)
    LUWA=np.where(MWU==1,4,LUWA)
    output_file=os.path.join(Output_dir,os.path.basename(WaPOR_LCC).replace('LCC','LUWA'))
    gis.CreateGeoTiff(output_file,LUWA,driver,NDV,xsize,ysize,GeoT,Projection)

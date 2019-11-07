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
import tempfile

def Rasterize_shapefile(fh, fh_ex,fh_out):
    """
    Rasterize a shapefile given an example tiff raster to get the pixelsizes and transformation info.
    
    Parameters
    ----------
    fh : str
        Filehandle to file to the shapefile to be rasterized.
    fh_eg : str
        Filehandle to file to the tiff raster example file.
    fh_out : str
        Filehandle for output.
        
    """

    shapefile = fh
    examplefile = fh_ex
    
    
    # Create temporary tif-file.
    temp_file = tempfile.NamedTemporaryFile(suffix='.tif').name
    
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(shapefile, 1)
    
    if not inDataSource:
        raise IOError("Could not open '%s'"%shapefile)
    
    inLayer = inDataSource.GetLayer()
    
    options = gdal.WarpOptions(cutlineDSName = shapefile,
                               cutlineLayer = inLayer.GetName(),
                               cropToCutline = False,
                               dstNodata = -9999,
                               )
    print(options)
    sourceds = gdal.Warp(temp_file, examplefile, options = options)
    if not sourceds:
        raise IOError("a problem in gdal.Warp '%s'"%shapefile)
    
    geot    = sourceds.GetGeoTransform()
    xsize   = sourceds.RasterXSize # columns
    ysize   = sourceds.RasterYSize # rows
    
    
    minX = geot[0]
    minY = geot[3] + ysize * geot[5]
    maxX = geot[0] + xsize * geot[1]
    maxY = geot[3]
    
    
    optionsProj = gdal.WarpOptions(
                               outputBounds = (minX, minY, maxX, maxY),
                               width = xsize,
                               height = ysize,
                               dstNodata = -9999,
                               options = ["GDALWARP_IGNORE_BAD_CUTLINE YES"],
                               )
    
    temp_fileP = tempfile.NamedTemporaryFile(suffix='.tif').name
    gdal.Warp(temp_fileP, temp_file, options = optionsProj)
    sourceds=gdal.Warp( fh_out,temp_fileP, options = optionsProj)
    sourceds=None
    os.remove(temp_file)
    os.remove(temp_fileP)

   
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

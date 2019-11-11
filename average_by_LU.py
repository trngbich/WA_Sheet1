# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 11:07:57 2019

@author: ntr002
"""

'''
input: xarray of 1 parameter, xarray of LU map (4 categories)

return: dataframe of yearly average values
'''
import xarray as xr
import pandas as pd


def Average_perLU(Var_xarray,LU_xarray):
    LU_dict={
            1:['PLU','Protected Landuse'],
            2:['ULU','Utilized Landuse'],
            3:['MLU','Modified Landuse'],
            4:['MWU','Managed Water Use'],
             }
    
    for i in LU_dict.keys():    
        V_LU=Var_xarray.where(LU_xarray==i)
        ts_V_LU=V_LU.mean(dim=['latitude','longitude']).to_dataframe()
        if i==1:
            ts_V=ts_V_LU.rename(columns={ts_V_LU.columns[0]:LU_dict[i][1]})   
        if i>1:
            ts_V_LU=ts_V_LU.rename(columns={ts_V_LU.columns[0]:LU_dict[i][1]})   
            ts_V=pd.merge(ts_V,ts_V_LU,left_index=True,right_index=True,how='outer')
    return ts_V

def Total_perLU(Var_xarray,LU_xarray):
    LU_dict={
            1:['PLU','Protected Landuse'],
            2:['ULU','Utilized Landuse'],
            3:['MLU','Modified Landuse'],
            4:['MWU','Managed Water Use'],
             }
    name=Var_xarray.name
    for i in LU_dict.keys():    
        V_LU=Var_xarray.where(LU_xarray==i)
        V_LU.name='{0}-{1}'.format(name,LU_dict[i][0])
        ts_V_LU=V_LU.sum(dim=['latitude','longitude']).to_dataframe()
        if i==1:
            ts_V=ts_V_LU.rename(columns={ts_V_LU.columns[0]:'-'.join([name,LU_dict[i][0]])})   
        if i>1:
            ts_V_LU=ts_V_LU.rename(columns={ts_V_LU.columns[0]:'-'.join([name,LU_dict[i][0]])})   
            ts_V=pd.merge(ts_V,ts_V_LU,left_index=True,right_index=True,how='outer')
    return ts_V
        


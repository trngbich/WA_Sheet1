# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 10:18:07 2018

@author: cmi001
"""
import os
import datetime
import calendar
import ogr
import numpy as np
from shapely.geometry import Polygon, Point
import shapefile


# from ogr cookbook https://pcjericks.github.io/py-gdalogr-cookbook/layers.html
def create_buffer(inputfn, output_bufferfn, buffer_dist):
    inputds = ogr.Open(inputfn)
    inputlyr = inputds.GetLayer()

    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(output_bufferfn):
#        os.remove(output_bufferfn)
        shpdriver.DeleteDataSource(output_bufferfn)
    outputBufferds = shpdriver.CreateDataSource(output_bufferfn)
    bufferlyr = outputBufferds.CreateLayer(output_bufferfn, geom_type=ogr.wkbPolygon)
    featureDefn = bufferlyr.GetLayerDefn()

    for feature in inputlyr:
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(buffer_dist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        bufferlyr.CreateFeature(outFeature)
        outFeature = None

def points_in_polygon(polyshp, pointcoords):
    polygon_r = shapefile.Reader(polyshp)
    polygon_shapes = polygon_r.shapes()
    shpfilePoints = []
    for shape in polygon_shapes:
        shpfilePoints = shape.points
    polygon_points = shpfilePoints
    polygon = Polygon(polygon_points)
    in_poly = []
    for coord in pointcoords:
        point = Point(coord)
        in_poly.append(polygon.contains(point))
    return np.where(np.array(in_poly))

### from bert
def convert_partial_year(number):
    year = int(number)
    d = datetime.timedelta(days=(number - year)*(365 + calendar.isleap(year)))
    day_one = datetime.date(year, 1, 1)
    date = d + day_one
    return date

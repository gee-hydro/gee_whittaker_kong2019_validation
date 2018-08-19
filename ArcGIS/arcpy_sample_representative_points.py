# -*- coding: utf-8 -*-
# Dongdong Kong, 20180520

import arcpy
from arcpy import env
from arcpy.sa import *
from arcpy import mapping
from os import path

mxd = mapping.MapDocument("CURRENT")
# layers = mapping.ListLayers(mxd)

indir = path.join(path.dirname(mxd.filePath), 'shp')
env.workspace = indir

# land cover in 0.05 deg
r = Raster("MCD12Q1_2010_006_20_full")  # no mask
r = Raster("MCD12Q1_2010_006_20_left")  # mask non-seasonality

number = '1000'  # how many random points generated?
shpfile = 'lc_poly_%s.shp' % (number)
shpfile_sv = 'lc_poly_%s_sv.shp' % (number)  # dissolved

file_st = 'st_' + number

print shpfile

# raster convert to polygon and simplify
arcpy.RasterToPolygon_conversion(r, shpfile, 'SIMPLIFY', 'VALUE')
arcpy.Dissolve_management(shpfile, shpfile_sv, 'GRIDCODE',
                          '#', 'MULTI_PART', 'DISSOLVE_LINES')

# rm water (IGBP == 0)
with arcpy.da.UpdateCursor(shpfile_sv, "GRIDCODE") as cursor:
    for row in cursor:
        if row[0] == 0.0:
            # print row
            cursor.deleteRow()

arcpy.RefreshActiveView()

# sample representative points
arcpy.CreateRandomPoints_management(
    indir, file_st, shpfile_sv, '0 0 250 250', number, '4 Kilometers', 'POINT', '0')

# https://pro.arcgis.com/zh-cn/pro-app/help/data/geodatabases/overview/arcgis-field-data-types.htm
arcpy.AddField_management(file_st, "site", "SHORT", 5)
arcpy.AddField_management(file_st, "IGBPcode", "SHORT", 2)

# add fields site and IGBPcode (CID + 1)
codeblock = (r'rec = 0 \n' +
             'def autoIncrement(): \n' +
             '    global rec \n' +
             '    pStart = 1 \n' +
             '    pInterval = 1 \n' +
             '    if (rec == 0): \n' +
             '        rec = pStart \n' +
             '    else: \n' +
             '        rec += pInterval \n' +
             '    return rec \n')
arcpy.CalculateField_management(file_st, 'site', 'autoIncrement()', 'PYTHON_9.3', codeblock)
arcpy.CalculateField_management(file_st, 'IGBPcode', '[CID]+1', 'VB', '#')

# delete unused fields
arcpy.DeleteField_management(file_st, ["CID"])


# lyr = mapping.ListLayers(mxd, file_st)[0]
# s = lyr.symbology
# s.classLabels = ['ENF', 'EBF', 'DNF', 'DBF', 'MF', 'CSH', 'OSH', 'WSA', 'SAV',
#                  'GRA', 'WET', 'CRO', 'URB', 'CNV', 'SNO', 'BSV']

# SAMPLE points in every kinds of IGBP to test phenofit
# for i in range(1, 18):
#     x = arcpy.sa.Con(r == i, r)
#     arcpy.Delete_management(shpfile, '#')

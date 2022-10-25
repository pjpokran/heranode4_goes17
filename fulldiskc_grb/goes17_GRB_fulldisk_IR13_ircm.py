#!/home/poker/miniconda3/envs/goes16_201710/bin/python

import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os,errno
#import os.rename
#import os.remove
import shutil
import sys
from PIL import Image

aoslogo = Image.open('/home/poker/uw-aoslogo.png')
aoslogoheight = aoslogo.size[1]
aoslogowidth = aoslogo.size[0]

# We need a float array between 0-1, rather than
# a uint8 array between 0-255
aoslogo = np.array(aoslogo).astype(np.float) / 255

band="13"
filechar=['AA','AB','AC','AD','AE','AF','AG','AH','AI','AJ','AK','AL','AM',
          'AN','AO','AP','AQ','AR','AS','AT','AU','AV','AW','AX','AY','AZ',
          'BA','BB','BC','BD','BE','BF','BG','BH','BI','BJ','BK','BL','BM',
          'BN','BO','BP','BQ','BR','BS','BT','BU','BV','BW','BX','BY','BZ']

#print(filechar[1])

prod_id = "TIRU"
#dt="201703051957"
dt = sys.argv[1]
#dt = 20172922011110

#f = netCDF4.Dataset("/weather/data/goes16/TIRC/"+band+"/"+dt+"_PAA.nc")
#f = netCDF4.Dataset("/weather/data/goes16/"+prod_id+"/"+band+"/"+dt+"_PAA.nc")

# ABI 2
#f = netCDF4.Dataset("/home/ldm/data/grb-test/OR_ABI-L1b-RadC-M3C02_G16_s20172991812207_e20172991814580_c20172991815015.nc")
# ABI 13
#f = netCDF4.Dataset("/home/ldm/data/grb-test/OR_ABI-L1b-RadF-M3C13_G16_s20172991845398_e20172991856176_c20172991856233.nc")
f = netCDF4.Dataset(dt)
#f = netCDF4.Dataset("/home/ldm/data/grb-test/OR_ABI-L1b-RadF-M3C02_G16_s20172921900306_e20172921911073_c20172921911108.nc")
#a = np.zeros(shape=(f.product_rows,f.product_columns))
#xa= np.zeros(shape=(f.product_columns))
#ya= np.zeros(shape=(f.product_rows))


print(f)

print("f.variables[Rad] ",f.variables["Rad"])

print("valid_pixel_count ",f.variables['valid_pixel_count'][0])
print("undersaturated_pixel_count ",f.variables['undersaturated_pixel_count'][0])
print("saturated_pixel_count ",f.variables['saturated_pixel_count'][0])
print("min_radiance_value_of_valid_pixels ",f.variables['min_radiance_value_of_valid_pixels'][0])
print("max_radiance_value_of_valid_pixels ",f.variables['max_radiance_value_of_valid_pixels'][0])
print("mean_radiance_value_of_valid_pixels ",f.variables['mean_radiance_value_of_valid_pixels'][0])
print("std_dev_radiance_value_of_valid_pixels ",f.variables['std_dev_radiance_value_of_valid_pixels'][0])
print("nominal_satellite_height ",f.variables['nominal_satellite_height'][0])
print("geospatial_lat_lon_extent ",f.variables['geospatial_lat_lon_extent'][0])
print("band_id ",f.variables['band_id'][0])
print("band_wavelength ",f.variables['band_wavelength'][0])
print("esun ",f.variables['esun'][0])
print("kappa0 ",f.variables['kappa0'][0])
print("earth_sun_distance_anomaly_in_AU ",f.variables['earth_sun_distance_anomaly_in_AU'][0])
print(" ")



# Convert radiance to Brightness Temp 
fk1 = f.variables['planck_fk1'][0]
fk2 = f.variables['planck_fk2'][0]
bc1 = f.variables['planck_bc1'][0]
bc2 = f.variables['planck_bc2'][0]
print("fk1 = ",fk1)
print("fk2 = ",fk2)
print("bc1 = ",bc1)
print("bc2 = ",bc2)
print(" ")
print("scale factor is ",f.variables['Rad'].scale_factor)
print("add_offset   is ",f.variables['Rad'].add_offset)
scalefactor = f.variables['Rad'].scale_factor
add_offset = f.variables['Rad'].add_offset
#data_var = f.variables['Rad'][0000:12000,0000:12000]
#data_var = f.variables['Rad'][:] *scalefactor + add_offset
data_var = f.variables['Rad'][:]
data_var = (fk2 / ( np.log((fk1 / data_var) + 1 )) - bc1) / bc2
print("BT values from y,x = 0-10, 0-10 (far northwest corner) ", data_var[0:10,0:10])
print("max of data_var ",np.max(data_var))
print("min of data_var ",np.min(data_var))
#data_var[data_var < 0] = 0.
#print("max of data_var ",np.max(data_var))
#print("min of data_var ",np.min(data_var))
a = data_var
print("max of data_var ",np.max(data_var))
print("min of data_var ",np.min(data_var))
#  T = [ fk2 / (alog((fk1 / Lλ) + 1)) - bc1 ] / bc2
print("radiance values from y,x = 0-10, 0-10 (far northwest corner) ", f.variables['Rad'][0:10,0:10])



print('x', f.variables['x'])
print('y', f.variables['y'])

# quit()


xa = f.variables['x'][:]
ya = f.variables['y'][:]

#    scale_factor: -1.4e-05
#    add_offset: 0.151865

print("xa",xa)
print("ya",ya)

#print("satellite height", f.variables['nominal_satellite_height'][0])
xa = xa*35785831
ya = ya*35785831
#ya = ya*f.variables['nominal_satellite_height'][0]/1000.

print("xa",xa)
print("ya",ya)

# quit()
# swap zeros for ones
# a[a==0.] = 1.

#print(np.average(a))
#if np.average(a) > 0.75:
#    quit()

import cartopy.crs as ccrs

##########################################################
## HACK for elipse and mask being drawn inside of the actual edge of the earth
#import shapely.geometry as sgeom
#def __init__(self, projection, satellite_height=35785831,
#             central_longitude=0.0, central_latitude=0.0,
#             false_easting=0, false_northing=0, globe=None):
#    proj4_params = [('proj', projection), ('lon_0', central_longitude),
#                    ('lat_0', central_latitude), ('h', satellite_height),
#                    ('x_0', false_easting), ('y_0', false_northing),
#                    ('units', 'm')]
#    super(ccrs._Satellite, self).__init__(proj4_params, globe=globe)
#
#    # TODO: Let the globe return the semimajor axis always.
##    a = np.float(self.globe.semimajor_axis or WGS84_SEMIMAJOR_AXIS )
#    a = np.float(self.globe.semimajor_axis or 6378137. )
#    b = np.float(self.globe.semiminor_axis or 6356752)
#    h = np.float(satellite_height)
#    max_x = h * np.arcsin(a / (a + h))
#    max_y = h * np.arcsin(b / (b + h))
#
#    coords = ccrs._ellipse_boundary(max_x, max_y,
#                                    false_easting, false_northing, 61)
#    self._boundary = sgeom.LinearRing(coords.T)
#    self._xlim = self._boundary.bounds[::2]
#    self._ylim = self._boundary.bounds[1::2]
#    self._threshold = np.diff(self._xlim)[0] * 0.02
#
#ccrs._Satellite.__init__ = __init__
#########################################################
##########################################################
## HACK for elipse and mask being drawn inside of the actual edge of the earth
#import math
#import shapely.geometry as sgeom
#def override_ellipse(self):
#    a = np.float(self.globe.semimajor_axis)
#    b = np.float(self.globe.semiminor_axis or a)
#    h = np.float(35785831.0)
#    max_x = 1.011 * h * math.atan(a / (a + h))
#    max_y = 1.011 * h * math.atan(b / (a + h))
#    coords = ccrs._ellipse_boundary(max_x, max_y, 0, 0, 61)
#    self._boundary = sgeom.LinearRing(coords.T)
#    self._xlim = self._boundary.bounds[::2]
#    self._ylim = self._boundary.bounds[1::2]
#    self._threshold = np.diff(self._xlim)[0] * 0.02
##########################################################

# Create a Globe specifying a spherical earth with the correct radius
#globe = ccrs.Globe(semimajor_axis=proj_var.semi_major,
#                   semiminor_axis=proj_var.semi_minor)
globe = ccrs.Globe(semimajor_axis=6378137.,semiminor_axis=6356752.)

#proj = ccrs.LambertConformal(central_longitude=proj_var.longitude_of_central_meridian,
#                             central_latitude=proj_var.latitude_of_projection_origin,
#                             standard_parallels=[proj_var.standard_parallel],
#                             globe=globe)

#proj = ccrs.Geostationary(central_longitude=-89.5, 
#                          satellite_height=35785831, 
proj = ccrs.Geostationary(central_longitude=f.variables['nominal_satellite_subpoint_lon'][0],
                          satellite_height=f.variables['nominal_satellite_height'][0] * 1000,
                          globe=globe,sweep_axis='x')

##########################################################
## HACK for elipse and mask being drawn inside of the actual edge of the earth
#override_ellipse(proj)
##########################################################

#f.variables['planck_bc1'][0]
print("y ",f.dimensions["y"].size)
print("x ",f.dimensions["x"].size)


image_rows=f.dimensions["y"].size
image_columns=f.dimensions["x"].size
namer_image_crop_top=0
namer_image_crop_bottom=-2400
namer_image_crop_left=400
namer_image_crop_right=-1400
#
namer_image_size_y=(image_rows+namer_image_crop_bottom-namer_image_crop_top)
namer_image_size_x=(image_columns+namer_image_crop_right-namer_image_crop_left)
#
print("namer image size")
print(namer_image_size_x, namer_image_size_y)
#
#wi_image_size_x=float(wi_image_size_x)/120.
#wi_image_size_y=float(wi_image_size_y)/120.
#namer_image_size_x=float(namer_image_size_x)/80.
#namer_image_size_y=float(namer_image_size_y)/80.
namer_image_size_x=15.1*1.325
namer_image_size_y=12.6*1.325


npac_image_crop_top=0
npac_image_crop_bottom=-3400
npac_image_crop_left=800
npac_image_crop_right=-1600
#
npac_image_size_y=(image_rows+npac_image_crop_bottom-npac_image_crop_top)
npac_image_size_x=(image_columns+npac_image_crop_right-npac_image_crop_left)
#
print("npac image size")
print(npac_image_size_x, npac_image_size_y)
#
#wi_image_size_x=float(wi_image_size_x)/120.
#wi_image_size_y=float(wi_image_size_y)/120.
#npac_image_size_x=float(npac_image_size_x)/80.
#npac_image_size_y=float(npac_image_size_y)/80.
npac_image_size_x=npac_image_size_x*.0075
npac_image_size_y=npac_image_size_y*.0075

ak_image_crop_top=10
ak_image_crop_bottom=-4800
ak_image_crop_left=1400
ak_image_crop_right=-2400
#
ak_image_size_y=(image_rows+ak_image_crop_bottom-ak_image_crop_top)
ak_image_size_x=(image_columns+ak_image_crop_right-ak_image_crop_left)
#
print("ak image size")
print(ak_image_size_x, ak_image_size_y)
#
#wi_image_size_x=float(wi_image_size_x)/120.
#wi_image_size_y=float(wi_image_size_y)/120.
#ak_image_size_x=float(ak_image_size_x)/80.
#ak_image_size_y=float(ak_image_size_y)/80.
ak_image_size_x=ak_image_size_x*.014
ak_image_size_y=ak_image_size_y*.014

#mw_image_crop_top=250
#mw_image_crop_bottom=-3200
#mw_image_crop_left=1750
#mw_image_crop_right=-2150
#
#mw_image_size_y=(image_rows+mw_image_crop_bottom-mw_image_crop_top)
#mw_image_size_x=(image_columns+mw_image_crop_right-mw_image_crop_left)
#
#print("mw image size")
#print(mw_image_size_x, mw_image_size_y)
#
#mw_image_size_x=float(mw_image_size_x)/150.
#mw_image_size_y=float(mw_image_size_y)/150.
#
#conus_image_crop_top=100
#conus_image_crop_bottom=-1800
#conus_image_crop_left=100
#conus_image_crop_right=-450
#
#conus_image_size_y=(image_rows+conus_image_crop_bottom-conus_image_crop_top)
#conus_image_size_x=(image_columns+conus_image_crop_right-conus_image_crop_left)
#
#print("conus image size")
#print(conus_image_size_x, conus_image_size_y)
#
#conus_image_size_x=float(conus_image_size_x)/300.
#conus_image_size_y=float(conus_image_size_y)/300.

# Create a new figure with size 10" by 10"
fig = plt.figure(figsize=(namer_image_size_x,namer_image_size_y),dpi=80.)
fig0 = plt.figure(figsize=(2.0,2.0))
#fig2 = plt.figure(figsize=(image_columns/160.,image_rows/160.),dpi=160.)
fig2 = plt.figure(figsize=(20.0,20.0),dpi=80.)    
fig20 = plt.figure(figsize=(2.0,2.0))
fig3 = plt.figure(figsize=(npac_image_size_x,npac_image_size_y),dpi=80.)
fig30 = plt.figure(figsize=(2.0,2.0))
fig4 = plt.figure(figsize=(ak_image_size_x,ak_image_size_y),dpi=80.)
fig40 = plt.figure(figsize=(2.0,2.0))
#fig2 = plt.figure(figsize=(image_columns/200.,image_rows/200.))
fig9 = plt.figure(figsize=(image_columns/78.,image_rows/78.))
#fig9 = plt.figure(figsize=(40,40))

# Put a single axes on this figure; set the projection for the axes to be our
# Lambert conformal projection
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.outline_patch.set_edgecolor('none')
ax0 = fig0.add_subplot(1, 1, 1, projection=proj)
ax0.outline_patch.set_edgecolor('none')
ax2 = fig2.add_subplot(1, 1, 1, projection=proj)
ax2.outline_patch.set_edgecolor('none')
ax20 = fig20.add_subplot(1, 1, 1, projection=proj)
ax20.outline_patch.set_edgecolor('none')
ax3 = fig3.add_subplot(1, 1, 1, projection=proj)
ax3.outline_patch.set_edgecolor('none')
ax30 = fig30.add_subplot(1, 1, 1, projection=proj)
ax30.outline_patch.set_edgecolor('none')
ax4 = fig4.add_subplot(1, 1, 1, projection=proj)
ax4.outline_patch.set_edgecolor('none')
ax40 = fig40.add_subplot(1, 1, 1, projection=proj)
ax40.outline_patch.set_edgecolor('none')
#ax3 = fig3.add_subplot(1, 1, 1, projection=proj)
ax9 = fig9.add_subplot(1, 1, 1, projection=proj)
ax9.outline_patch.set_edgecolor('none')


import matplotlib as mpl

cdict = {'red': ((0.0, 0.0, 0.0),
                 (.001, 1.00, 1.00),
                 (.107, 1.00, 1.00),
                 (.113, 0.498, 0.498),
                 (.173, 1.00, 1.00),
                 (.179, 0.902, 0.902),
                 (.227, 0.102, 0.102),
                 (.233, 0.00, 0.00),
                 (.287, 0.902, 0.902),
                 (.293, 1.00, 1.00),
                 (.346, 1.00, 1.00),
                 (.352, 1.00, 1.00),
                 (.406, 0.101, 0.101),
                (.412, 0.00, 0.00),
                 (.481, 0.00, 0.00),
                 (.484, 0.00, 0.00),
                 (.543, 0.00, 0.00),
                 (.546, 0.773, 0.773),
                 (.994, 0.012, 0.012),
                 (.997, 0.004, 0.004),
                 (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.0, 0.0),
                 (.001, 1.00, 1.00),
                 (.107, 1.00, 1.00),
                 (.113, 0.00, 0.00),
                 (.173, 0.498, 0.498),
                 (.179, 0.902, 0.902),
                 (.227, 0.102, 0.102),
                 (.233, 0.00, 0.00),
                 (.287, 0.00, 0.00),
                 (.293, 0.00, 0.00),
                 (.346, 0.902, 0.902),
                 (.352, 1.00, 1.00),
                 (.406, 1.00, 1.00),
                 (.412, 1.00, 1.00),
                 (.481, 0.00, 0.00),
                 (.484, 0.00, 0.00),
                 (.543, 1.00, 1.00),
                 (.546, 0.773, 0.773),
                 (.994, 0.012, 0.012),
                 (.997, 0.004, 0.004),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.00, 0.00),
                 (.001, 1.00, 1.00),
                 (.107, 0.00, 0.00),
                 (.113, 0.498, 0.498),
                 (.173, 0.786, 0.786),
                 (.179, 0.902, 0.902),
                 (.227, 0.102, 0.102),
                 (.233, 0.00, 0.00),
                 (.287, 0.00, 0.00),
                 (.293, 0.00, 0.00),
                 (.346, 0.00, 0.00),
                 (.352, 0.00, 0.00),
                 (.406, 0.00, 0.00),
                 (.412, 0.00, 0.00),
                 (.481, 0.451, 0.451),
                 (.484, 0.451, 0.451),
                 (.543, 1.00, 1.00),
                 (.546, 0.773, 0.773),
                 (.994, 0.012, 0.012),
                 (.997, 0.004, 0.004),
                  (1.0, 0.0, 0.0))}


my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,2048)

# Plot the data using a simple greyscale colormap (with black for low values);
# set the colormap to extend over a range of values from 140 to 255.
# Note, we save the image returned by imshow for later...
#im = ax.imshow(data_var[:], extent=(x[0], x[-1], y[0], y[-1]), origin='upper',
#               cmap='Greys_r', norm=plt.Normalize(0, 256))
#im = ax.imshow(data_var[:], extent=(x[0], x[-1], y[0], y[-1]), origin='upper',
#im = ax.imshow(a[:], extent=(xa[0], xa[-1], ya[-1], ya[0]), origin='upper',
#               cmap='Greys_r')
#im = ax.imshow(a[250:-3000,2000:-2000], extent=(xa[2000],xa[-2000],ya[-3000],ya[250]), origin='upper',cmap='Greys_r')
#im = ax.imshow(a[250:-3500,2500:-2200], extent=(xa[2500],xa[-2200],ya[-3500],ya[250]), origin='upper',cmap='Greys_r')
#im = ax.imshow(data_var[:], extent=(x[0], x[-1], y[0], y[-1]), origin='upper')
#im = ax2.imshow(a[:], extent=(xa[1],xa[-1],ya[-1],ya[1]), origin='upper', cmap='Greys_r')

im = ax.imshow(a[namer_image_crop_top:namer_image_crop_bottom,namer_image_crop_left:namer_image_crop_right], extent=(xa[namer_image_crop_left],xa[namer_image_crop_right],ya[namer_image_crop_bottom],ya[namer_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
im0 = ax0.imshow(a[namer_image_crop_top:namer_image_crop_bottom,namer_image_crop_left:namer_image_crop_right], extent=(xa[namer_image_crop_left],xa[namer_image_crop_right],ya[namer_image_crop_bottom],ya[namer_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
#im = ax2.imshow(a[mw_image_crop_top:mw_image_crop_bottom,mw_image_crop_left:mw_image_crop_right], extent=(xa[mw_image_crop_left],xa[mw_image_crop_right],ya[mw_image_crop_bottom],ya[mw_image_crop_top]), origin='upper',cmap='Greys_r', vmin=0., vmax=1.)
im2 = ax2.imshow(a[:], extent=(xa[0],xa[-1],ya[-1],ya[0]), origin='upper', cmap=my_cmap, vmin=162., vmax=330.)
im20 = ax20.imshow(a[:], extent=(xa[0],xa[-1],ya[-1],ya[0]), origin='upper', cmap=my_cmap, vmin=162., vmax=330.)
im3 = ax3.imshow(a[npac_image_crop_top:npac_image_crop_bottom,npac_image_crop_left:npac_image_crop_right], extent=(xa[npac_image_crop_left],xa[npac_image_crop_right],ya[npac_image_crop_bottom],ya[npac_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
im30 = ax30.imshow(a[npac_image_crop_top:npac_image_crop_bottom,npac_image_crop_left:npac_image_crop_right], extent=(xa[npac_image_crop_left],xa[npac_image_crop_right],ya[npac_image_crop_bottom],ya[npac_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
im4 = ax4.imshow(a[ak_image_crop_top:ak_image_crop_bottom,ak_image_crop_left:ak_image_crop_right], extent=(xa[ak_image_crop_left],xa[ak_image_crop_right],ya[ak_image_crop_bottom],ya[ak_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
im40 = ax40.imshow(a[ak_image_crop_top:ak_image_crop_bottom,ak_image_crop_left:ak_image_crop_right], extent=(xa[ak_image_crop_left],xa[ak_image_crop_right],ya[ak_image_crop_bottom],ya[ak_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
#im = ax3.imshow(a[conus_image_crop_top:conus_image_crop_bottom,conus_image_crop_left:conus_image_crop_right], extent=(xa[conus_image_crop_left],xa[conus_image_crop_right],ya[conus_image_crop_bottom],ya[conus_image_crop_top]), origin='upper',cmap='Greys_r', vmin=0., vmax=1.)
im9 = ax9.imshow(a[:], extent=(xa[0],xa[-1],ya[-1],ya[0]), origin='upper', cmap=my_cmap, vmin=162., vmax=330.)

ax.coastlines(resolution='50m', color='green')
ax0.coastlines(resolution='50m', color='green')
ax2.coastlines(resolution='50m', color='green')
ax20.coastlines(resolution='50m', color='green')
ax3.coastlines(resolution='50m', color='green')
ax30.coastlines(resolution='50m', color='green')
ax4.coastlines(resolution='50m', color='green')
ax40.coastlines(resolution='50m', color='green')
ax9.coastlines(resolution='50m', color='green')

import cartopy.feature as cfeat

# Add country borders with a thick line.
ax.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax0.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax2.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax20.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax3.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax30.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax4.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax40.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
#ax3.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
ax9.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')

# Set up a feature for the state/province lines. Tell cartopy not to fill in the polygons
state_boundaries = cfeat.NaturalEarthFeature(category='cultural',
                                             name='admin_1_states_provinces_lakes',
                                             scale='50m', facecolor='none', edgecolor='red')

# Add the feature with dotted lines, denoted by ':'
ax.add_feature(state_boundaries, linestyle=':')
ax0.add_feature(state_boundaries, linestyle=':')
ax2.add_feature(state_boundaries, linestyle=':')
ax20.add_feature(state_boundaries, linestyle=':')
ax3.add_feature(state_boundaries, linestyle=':')
ax30.add_feature(state_boundaries, linestyle=':')
ax4.add_feature(state_boundaries, linestyle=':')
ax40.add_feature(state_boundaries, linestyle=':')
#ax3.add_feature(state_boundaries, linestyle=':')
ax9.add_feature(state_boundaries, linestyle=':')

# axes for namer
cbaxes1 = fig.add_axes([0.135,0.14,0.755,0.02])
cbar1 = fig.colorbar(im, cax=cbaxes1, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar1.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar1.ax.xaxis.set_ticks_position('top')
cbar1.ax.xaxis.set_label_position('top')

# axes for fulldisk
#cbaxes2 = fig2.add_axes([0.135,0.12,0.755,0.02])
cbaxes2 = fig2.add_axes([0.135,0.12,0.255,0.005])
cbar2 = fig2.colorbar(im2, cax=cbaxes2, orientation='horizontal')
font_size = 10
#cbar2.set_label('Brightness Temperature (K)',size=18)
cbar2.ax.tick_params(labelsize=font_size, labelcolor='black')
cbar2.ax.xaxis.set_ticks_position('top')
cbar2.ax.xaxis.set_label_position('top')

# axes for npac
cbaxes3 = fig3.add_axes([0.135,0.14,0.755,0.02])
cbar3 = fig3.colorbar(im, cax=cbaxes3, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar3.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar3.ax.xaxis.set_ticks_position('top')
cbar3.ax.xaxis.set_label_position('top')

# axes for ak
cbaxes4 = fig4.add_axes([0.135,0.14,0.755,0.02])
cbar4 = fig4.colorbar(im, cax=cbaxes4, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar4.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar4.ax.xaxis.set_ticks_position('top')
cbar4.ax.xaxis.set_label_position('top')

# axes for full
cbaxes9 = fig9.add_axes([0.135,0.15,0.755,0.02])
cbar9 = fig9.colorbar(im9, cax=cbaxes9, orientation='horizontal')
font_size = 18
#cbar9.set_label('Brightness Temperature (K)',size=18)
cbar9.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar9.ax.xaxis.set_ticks_position('top')
cbar9.ax.xaxis.set_label_position('top')




# Redisplay modified figure
#fig
#fig2

import datetime

#time_var = f.date_created
time_var = f.time_coverage_start

#jyr = time_var[0:4]
#jday = time_var[4:7]
#print(jday)
iyear = time_var[0:4]
print("iyear ",iyear)
imonth = time_var[5:7]
print("imonth ",imonth)
import calendar
cmonth = calendar.month_abbr[int(imonth)]
print("cmonth ",cmonth)
iday = time_var[8:10]
print("iday ",iday)
itime = time_var[11:19]
itimehr = time_var[11:13]
itimemn = time_var[14:16]

#date = datetime.datetime(int(jyr), 1, 1) + datetime.timedelta(int(jday)-1)

ctime_string = iyear +' '+cmonth+' '+iday+'  '+itime+' GMT'
ctime_file_string = iyear + imonth + iday + itimehr + itimemn

#time_string = 'GOES-16 Band 14 (LW IR) valid %s'%time_var
#time_string = 'GOES-16 Band 14 (LW IR) valid %s %s %s %s'%iyear %cmonth %iday %itime
time_string = 'GOES-17 Band 13\nClean LW IR Window\n%s'%ctime_string
print(time_string)

from matplotlib import patheffects
outline_effect = [patheffects.withStroke(linewidth=2, foreground='black')]

#2017/065 20:04:00:30
text = ax.text(0.005, 0.89, time_string,
    horizontalalignment='left', transform = ax.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text.set_path_effects(outline_effect)
#
text2 = ax2.text(0.005, 0.92, time_string,
    horizontalalignment='left', transform = ax2.transAxes,
    color='yellow', fontsize='14', weight='bold')

text2.set_path_effects(outline_effect)

text3 = ax3.text(0.005, 0.89, time_string,
    horizontalalignment='left', transform = ax3.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text3.set_path_effects(outline_effect)
#
text4 = ax4.text(0.005, 0.80, time_string,
    horizontalalignment='left', transform = ax4.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text4.set_path_effects(outline_effect)
#
#text = ax3.text(0.50, 0.90, time_string,
#    horizontalalignment='center', transform = ax3.transAxes,
#    color='yellow', fontsize='large', weight='bold')

text9 = ax9.text(0.50, 0.97, time_string,
    horizontalalignment='center', transform = ax9.transAxes,
    color='yellow', fontsize='large', weight='bold')

text9.set_path_effects(outline_effect)



filename1="/whirlwind/goes17/irc13m/namer/"+ctime_file_string+"_namer.jpg"
filename0="/whirlwind/goes17/irc13m/namer/latest_namer_1s.jpg"
filename2="/whirlwind/goes17/irc13m/fulldisk/"+ctime_file_string+"_fulldisk.jpg"
filename20="/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_1s.jpg"
filename3="/whirlwind/goes17/irc13m/npac/"+ctime_file_string+"_npac.jpg"
filename30="/whirlwind/goes17/irc13m/npac/latest_npac_1s.jpg"
filename4="/whirlwind/goes17/irc13m/ak/"+ctime_file_string+"_ak.jpg"
filename40="/whirlwind/goes17/irc13m/ak/latest_ak_1s.jpg"
#filename3="/whirlwind/goes17/vis_sqrt/conus/"+dt+"_conus.jpg"
filename9="/whirlwind/goes17/irc13m/fulldisk_full/"+ctime_file_string+"_fulldisk_full.jpg"

fig.figimage(aoslogo,  0, fig.bbox.ymax*0.77 - aoslogoheight    , zorder=10)
fig2.figimage(aoslogo,  0, fig2.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig3.figimage(aoslogo,  0, fig3.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig4.figimage(aoslogo,  0, fig4.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig9.figimage(aoslogo,  0, fig9.bbox.ymax*0.77 - aoslogoheight  , zorder=10)

fig.savefig(filename1, bbox_inches='tight', pad_inches=0)
fig0.savefig(filename0, bbox_inches='tight', pad_inches=0)
fig2.savefig(filename2, bbox_inches='tight', pad_inches=0)
fig20.savefig(filename20, bbox_inches='tight', pad_inches=0)
fig3.savefig(filename3, bbox_inches='tight', pad_inches=0)
fig30.savefig(filename30, bbox_inches='tight', pad_inches=0)
fig4.savefig(filename4, bbox_inches='tight', pad_inches=0)
fig40.savefig(filename40, bbox_inches='tight', pad_inches=0)
#fig2.savefig(filename2jpg, bbox_inches='tight')
#fig3.savefig(filename3, bbox_inches='tight')
fig9.savefig(filename9, bbox_inches='tight', pad_inches=0)

#quit()

#import os.rename    # os.rename(src,dest)
#import os.remove    # os.remove path
#import shutil.copy  # shutil.copy(src, dest)


def silentremove(filename):
    try: 
        os.remove(filename)
    except OSError:
        pass
    
def silentrename(filename1, filename2):
    try: 
        os.rename(filename1, filename2)
    except OSError:
        pass
    
silentremove("/whirlwind/goes17/irc13m/namer/latest_namer_24.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_23.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_24.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_22.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_23.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_21.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_22.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_20.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_21.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_19.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_20.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_18.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_19.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_17.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_18.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_16.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_17.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_15.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_16.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_14.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_15.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_13.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_14.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_12.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_13.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_11.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_12.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_10.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_11.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_9.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_10.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_8.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_9.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_7.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_8.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_6.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_7.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_5.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_6.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_4.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_5.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_3.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_4.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_2.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_3.jpg")
silentrename("/whirlwind/goes17/irc13m/namer/latest_namer_1.jpg", "/whirlwind/goes17/irc13m/namer/latest_namer_2.jpg")
    
shutil.copy(filename1, "/whirlwind/goes17/irc13m/namer/latest_namer_1.jpg")


silentremove("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_24.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_23.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_24.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_22.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_23.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_21.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_22.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_20.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_21.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_19.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_20.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_18.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_19.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_17.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_18.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_16.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_17.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_15.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_16.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_14.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_15.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_13.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_14.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_12.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_13.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_11.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_12.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_10.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_11.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_9.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_10.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_8.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_9.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_7.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_8.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_6.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_7.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_5.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_6.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_4.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_5.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_3.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_4.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_2.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_3.jpg")
silentrename("/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_1.jpg", "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_2.jpg")
    
shutil.copy(filename2, "/whirlwind/goes17/irc13m/fulldisk/latest_fulldisk_1.jpg")

silentremove("/whirlwind/goes17/irc13m/npac/latest_npac_24.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_23.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_24.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_22.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_23.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_21.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_22.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_20.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_21.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_19.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_20.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_18.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_19.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_17.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_18.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_16.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_17.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_15.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_16.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_14.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_15.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_13.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_14.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_12.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_13.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_11.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_12.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_10.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_11.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_9.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_10.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_8.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_9.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_7.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_8.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_6.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_7.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_5.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_6.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_4.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_5.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_3.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_4.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_2.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_3.jpg")
silentrename("/whirlwind/goes17/irc13m/npac/latest_npac_1.jpg", "/whirlwind/goes17/irc13m/npac/latest_npac_2.jpg")
    
shutil.copy(filename3, "/whirlwind/goes17/irc13m/npac/latest_npac_1.jpg")

silentremove("/whirlwind/goes17/irc13m/ak/latest_ak_24.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_23.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_24.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_22.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_23.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_21.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_22.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_20.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_21.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_19.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_20.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_18.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_19.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_17.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_18.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_16.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_17.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_15.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_16.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_14.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_15.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_13.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_14.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_12.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_13.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_11.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_12.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_10.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_11.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_9.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_10.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_8.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_9.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_7.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_8.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_6.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_7.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_5.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_6.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_4.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_5.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_3.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_4.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_2.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_3.jpg")
silentrename("/whirlwind/goes17/irc13m/ak/latest_ak_1.jpg", "/whirlwind/goes17/irc13m/ak/latest_ak_2.jpg")
    
shutil.copy(filename4, "/whirlwind/goes17/irc13m/ak/latest_ak_1.jpg")

shutil.copy(filename9, "/whirlwind/goes17/irc13m/fulldisk_full/latest_fulldisk_full_1.jpg")

quit()

os.remove("/whirlwind/goes16/vis/conus/latest_conus_72.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_71.jpg", "/whirlwind/goes16/vis/conus/latest_conus_72.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_70.jpg", "/whirlwind/goes16/vis/conus/latest_conus_71.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_69.jpg", "/whirlwind/goes16/vis/conus/latest_conus_70.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_68.jpg", "/whirlwind/goes16/vis/conus/latest_conus_69.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_67.jpg", "/whirlwind/goes16/vis/conus/latest_conus_68.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_66.jpg", "/whirlwind/goes16/vis/conus/latest_conus_67.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_65.jpg", "/whirlwind/goes16/vis/conus/latest_conus_66.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_64.jpg", "/whirlwind/goes16/vis/conus/latest_conus_65.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_63.jpg", "/whirlwind/goes16/vis/conus/latest_conus_64.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_62.jpg", "/whirlwind/goes16/vis/conus/latest_conus_63.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_61.jpg", "/whirlwind/goes16/vis/conus/latest_conus_62.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_60.jpg", "/whirlwind/goes16/vis/conus/latest_conus_61.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_59.jpg", "/whirlwind/goes16/vis/conus/latest_conus_60.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_58.jpg", "/whirlwind/goes16/vis/conus/latest_conus_59.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_57.jpg", "/whirlwind/goes16/vis/conus/latest_conus_58.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_56.jpg", "/whirlwind/goes16/vis/conus/latest_conus_57.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_55.jpg", "/whirlwind/goes16/vis/conus/latest_conus_56.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_54.jpg", "/whirlwind/goes16/vis/conus/latest_conus_55.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_53.jpg", "/whirlwind/goes16/vis/conus/latest_conus_54.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_52.jpg", "/whirlwind/goes16/vis/conus/latest_conus_53.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_51.jpg", "/whirlwind/goes16/vis/conus/latest_conus_52.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_50.jpg", "/whirlwind/goes16/vis/conus/latest_conus_51.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_49.jpg", "/whirlwind/goes16/vis/conus/latest_conus_50.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_48.jpg", "/whirlwind/goes16/vis/conus/latest_conus_49.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_47.jpg", "/whirlwind/goes16/vis/conus/latest_conus_48.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_46.jpg", "/whirlwind/goes16/vis/conus/latest_conus_47.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_45.jpg", "/whirlwind/goes16/vis/conus/latest_conus_46.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_44.jpg", "/whirlwind/goes16/vis/conus/latest_conus_45.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_43.jpg", "/whirlwind/goes16/vis/conus/latest_conus_44.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_42.jpg", "/whirlwind/goes16/vis/conus/latest_conus_43.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_41.jpg", "/whirlwind/goes16/vis/conus/latest_conus_42.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_40.jpg", "/whirlwind/goes16/vis/conus/latest_conus_41.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_39.jpg", "/whirlwind/goes16/vis/conus/latest_conus_40.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_38.jpg", "/whirlwind/goes16/vis/conus/latest_conus_39.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_37.jpg", "/whirlwind/goes16/vis/conus/latest_conus_38.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_36.jpg", "/whirlwind/goes16/vis/conus/latest_conus_37.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_35.jpg", "/whirlwind/goes16/vis/conus/latest_conus_36.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_34.jpg", "/whirlwind/goes16/vis/conus/latest_conus_35.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_33.jpg", "/whirlwind/goes16/vis/conus/latest_conus_34.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_32.jpg", "/whirlwind/goes16/vis/conus/latest_conus_33.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_31.jpg", "/whirlwind/goes16/vis/conus/latest_conus_32.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_30.jpg", "/whirlwind/goes16/vis/conus/latest_conus_31.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_29.jpg", "/whirlwind/goes16/vis/conus/latest_conus_30.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_28.jpg", "/whirlwind/goes16/vis/conus/latest_conus_29.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_27.jpg", "/whirlwind/goes16/vis/conus/latest_conus_28.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_26.jpg", "/whirlwind/goes16/vis/conus/latest_conus_27.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_25.jpg", "/whirlwind/goes16/vis/conus/latest_conus_26.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_24.jpg", "/whirlwind/goes16/vis/conus/latest_conus_25.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_23.jpg", "/whirlwind/goes16/vis/conus/latest_conus_24.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_22.jpg", "/whirlwind/goes16/vis/conus/latest_conus_23.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_21.jpg", "/whirlwind/goes16/vis/conus/latest_conus_22.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_20.jpg", "/whirlwind/goes16/vis/conus/latest_conus_21.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_19.jpg", "/whirlwind/goes16/vis/conus/latest_conus_20.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_18.jpg", "/whirlwind/goes16/vis/conus/latest_conus_19.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_17.jpg", "/whirlwind/goes16/vis/conus/latest_conus_18.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_16.jpg", "/whirlwind/goes16/vis/conus/latest_conus_17.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_15.jpg", "/whirlwind/goes16/vis/conus/latest_conus_16.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_14.jpg", "/whirlwind/goes16/vis/conus/latest_conus_15.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_13.jpg", "/whirlwind/goes16/vis/conus/latest_conus_14.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_12.jpg", "/whirlwind/goes16/vis/conus/latest_conus_13.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_11.jpg", "/whirlwind/goes16/vis/conus/latest_conus_12.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_10.jpg", "/whirlwind/goes16/vis/conus/latest_conus_11.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_9.jpg", "/whirlwind/goes16/vis/conus/latest_conus_10.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_8.jpg", "/whirlwind/goes16/vis/conus/latest_conus_9.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_7.jpg", "/whirlwind/goes16/vis/conus/latest_conus_8.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_6.jpg", "/whirlwind/goes16/vis/conus/latest_conus_7.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_5.jpg", "/whirlwind/goes16/vis/conus/latest_conus_6.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_4.jpg", "/whirlwind/goes16/vis/conus/latest_conus_5.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_3.jpg", "/whirlwind/goes16/vis/conus/latest_conus_4.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_2.jpg", "/whirlwind/goes16/vis/conus/latest_conus_3.jpg")
os.rename("/whirlwind/goes16/vis/conus/latest_conus_1.jpg", "/whirlwind/goes16/vis/conus/latest_conus_2.jpg")

shutil.copy(filename3, "/whirlwind/goes16/vis/conus/latest_conus_1.jpg")


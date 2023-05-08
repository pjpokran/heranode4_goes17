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
import glob
from PIL import Image

aoslogo = Image.open('/home/poker/uw-aoslogo.png')
aoslogoheight = aoslogo.size[1]
aoslogowidth = aoslogo.size[0]

# We need a float array between 0-1, rather than
# a uint8 array between 0-255
aoslogo = np.array(aoslogo).astype(np.float) / 255

band="09"
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
# ABI 09
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

# Create a Globe specifying a spherical earth with the correct radius
globe = ccrs.Globe(semimajor_axis=6378137.,semiminor_axis=6356752.)

proj = ccrs.Geostationary(central_longitude=f.variables['nominal_satellite_subpoint_lon'][0],
                          satellite_height=f.variables['nominal_satellite_height'][0] * 1000,
                          globe=globe,sweep_axis='x')

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
ak_image_size_x=ak_image_size_x*.014
ak_image_size_y=ak_image_size_y*.014

import matplotlib as mpl
import cartopy.feature as cfeat

# Set up a feature for the state/province lines. Tell cartopy not to fill in the polygons
state_boundaries = cfeat.NaturalEarthFeature(category='cultural',
                                             name='admin_1_states_provinces_lakes',
                                             scale='50m', facecolor='none', edgecolor='red')

# Set up a feature for the lat/lon lines at 15 deg. Tell cartopy not to fill in the polgons
lat_lon_lines_15deg = cfeat.NaturalEarthFeature(category='physical',
                                             name='graticules_15',
                                             scale='50m', facecolor='none', edgecolor='magenta')

from matplotlib import patheffects
outline_effect = [patheffects.withStroke(linewidth=2, foreground='black')]


# Create color table for color enhanced version
cdict = {'red': ((0.0, 0.0, 0.0),
                 (0.290, 0.263, .263),
                 (0.385, 1.0, 1.0),
                 (0.475, 0.443, .443),
                 (0.515, 0.0, 0.0),
                 (0.575, 1.0, 1.0),
                 (0.664, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.0, 0.0),
                   (0.290, .513, .513),
                   (0.385, 1.0, 1.0),
                   (0.475, .443, .443),
                   (0.515, 0., 0.0),
                   (0.575, 1.0, 1.0),
                   (0.664, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.0, 0.0),
                  (0.290, .137, .137),
                  (0.385, 1.0, 1.0),
                  (0.475,0.694, 0.694),
                  (0.515, .451, .451),
                  (0.552, 0.0, 0.0),
                  (0.664, 0.0, 0.0),
                  (1.0, 0.0, 0.0))}

my_cmap = mpl.colors.LinearSegmentedColormap('my_colormap',cdict,2048)

# Color table for b/w enhanced
cdict2 = {'red': ((0.0, 1.0, 1.0),
                 (0.29, 1.00, 1.00),
                 (0.61, 0.0, 0.0),
                 (1.0, 0.0, 0.0)),
         'green': ((0.0, 1.0, 1.0),
                 (0.29, 1.00, 1.00),
                 (0.61, 0.0, 0.0),
                 (1.0, 0.0, 0.0)),
         'blue': ((0.0, 1.0, 1.0),
                 (0.29, 1.00, 1.00),
                 (0.61, 0.0, 0.0),
                 (1.0, 0.0, 0.0))}
my_cmap2 = mpl.colors.LinearSegmentedColormap('my_colormap2',cdict2,2048)


import datetime

time_var = f.time_coverage_start

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


ctime_string = iyear +' '+cmonth+' '+iday+'  '+itime+' GMT'
ctime_file_string = iyear + imonth + iday + itimehr + itimemn

if f.platform_ID == "G17":
    time_string = 'GOES-17 Band 9\nMid-level Water Vapor\n%s'%ctime_string
elif f.platform_ID == "G18":
    time_string = 'GOES-18 Band 9\nMid-level Water Vapor\n%s'%ctime_string
else:
    time_string = 'GOES-West Band 9\nMid-level Water Vapor\n%s'%ctime_string
print(time_string)

# START NAMER
# Fig for color namer
fig = plt.figure(figsize=(namer_image_size_x,namer_image_size_y),dpi=80.)
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.outline_patch.set_edgecolor('none')
im = ax.imshow(a[namer_image_crop_top:namer_image_crop_bottom,namer_image_crop_left:namer_image_crop_right], extent=(xa[namer_image_crop_left],xa[namer_image_crop_right],ya[namer_image_crop_bottom],ya[namer_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
ax.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax.add_feature(state_boundaries, linestyle=':')
ax.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for namer
cbaxes1 = fig.add_axes([0.135,0.14,0.755,0.02])
cbar1 = fig.colorbar(im, cax=cbaxes1, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar1.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar1.ax.xaxis.set_ticks_position('top')
cbar1.ax.xaxis.set_label_position('top')
text = ax.text(0.005, 0.89, time_string,
    horizontalalignment='left', transform = ax.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text.set_path_effects(outline_effect)
#
filename1="/whirlwind/goes17/wvc/namer/"+ctime_file_string+"_namer.jpg"
#filename1="/dev/shm/wvc_"+ctime_file_string+"_namer.jpg"
fig.figimage(aoslogo,  0, fig.bbox.ymax*0.77 - aoslogoheight    , zorder=10)
fig.savefig(filename1, bbox_inches='tight', pad_inches=0)
shutil.copy(filename1, "/whirlwind/goes17/wvc/namer/latest_namer_1.jpg")
plt.close(fig)

# NAMER
file_list = glob.glob('/whirlwind/goes17/wvc/namer/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 1 week // every 3rd image (337 by 3)
thefile = open('/whirlwind/goes17/wvc/namer_week_temp.list', 'w')
thelist = file_list[-1009::3]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/namer_week_temp.list','/whirlwind/goes17/wvc/goes17_namer_loop_week.list')

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wvc/namer_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/namer_24_temp.list','/whirlwind/goes17/wvc/goes17_namer_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wvc/namer_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/namer_60_temp.list','/whirlwind/goes17/wvc/goes17_namer_loop_60.list')


# Fig for b/w namer
fig = plt.figure(figsize=(namer_image_size_x,namer_image_size_y),dpi=80.)
ax = fig.add_subplot(1, 1, 1, projection=proj)
ax.outline_patch.set_edgecolor('none')
im = ax.imshow(a[namer_image_crop_top:namer_image_crop_bottom,namer_image_crop_left:namer_image_crop_right], extent=(xa[namer_image_crop_left],xa[namer_image_crop_right],ya[namer_image_crop_bottom],ya[namer_image_crop_top]), origin='upper',cmap=my_cmap2, vmin=162., vmax=330.)
ax.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax.add_feature(state_boundaries, linestyle=':')
ax.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for namer
cbaxes1 = fig.add_axes([0.135,0.14,0.755,0.02])
cbar1 = fig.colorbar(im, cax=cbaxes1, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar1.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar1.ax.xaxis.set_ticks_position('top')
cbar1.ax.xaxis.set_label_position('top')
text = ax.text(0.005, 0.89, time_string,
    horizontalalignment='left', transform = ax.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text.set_path_effects(outline_effect)
#
filename1="/whirlwind/goes17/wv/namer/"+ctime_file_string+"_namer.jpg"
#filename1="/dev/shm/wv_"+ctime_file_string+"_namer.jpg"
fig.figimage(aoslogo,  0, fig.bbox.ymax*0.77 - aoslogoheight    , zorder=10)
fig.savefig(filename1, bbox_inches='tight', pad_inches=0)
shutil.copy(filename1, "/whirlwind/goes17/wv/namer/latest_namer_1.jpg")
plt.close(fig)

# NAMER
file_list = glob.glob('/whirlwind/goes17/wv/namer/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wv/namer_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/namer_24_temp.list','/whirlwind/goes17/wv/goes17_namer_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wv/namer_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/namer_60_temp.list','/whirlwind/goes17/wv/goes17_namer_loop_60.list')



# END NAMER

# START Fulldisk
# Fig for color fulldisk
fig2 = plt.figure(figsize=(20.0,20.0),dpi=80.)    
ax2 = fig2.add_subplot(1, 1, 1, projection=proj)
ax2.outline_patch.set_edgecolor('none')
im2 = ax2.imshow(a[:], extent=(xa[0],xa[-1],ya[-1],ya[0]), origin='upper', cmap=my_cmap, vmin=162., vmax=330.)
ax2.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax2.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax2.add_feature(state_boundaries, linestyle=':')
ax2.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for fulldisk
#cbaxes2 = fig2.add_axes([0.135,0.12,0.755,0.02])
cbaxes2 = fig2.add_axes([0.135,0.12,0.255,0.005])
cbar2 = fig2.colorbar(im2, cax=cbaxes2, orientation='horizontal')
font_size = 10
#cbar2.set_label('Brightness Temperature (K)',size=18)
cbar2.ax.tick_params(labelsize=font_size, labelcolor='black')
cbar2.ax.xaxis.set_ticks_position('top')
cbar2.ax.xaxis.set_label_position('top')
text2 = ax2.text(0.005, 0.92, time_string,
    horizontalalignment='left', transform = ax2.transAxes,
    color='yellow', fontsize='14', weight='bold')

text2.set_path_effects(outline_effect)

filename2="/whirlwind/goes17/wvc/fulldisk/"+ctime_file_string+"_fulldisk.jpg"
#filename2="/dev/shm/wvc_"+ctime_file_string+"_fulldisk.jpg"
fig2.figimage(aoslogo,  0, fig2.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig2.savefig(filename2, bbox_inches='tight', pad_inches=0)
shutil.copy(filename2, "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_1.jpg")
plt.close(fig2)

# FULLDISK Color
file_list = glob.glob('/whirlwind/goes17/wvc/fulldisk/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wvc/fulldisk_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/fulldisk_24_temp.list','/whirlwind/goes17/wvc/goes17_fulldisk_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wvc/fulldisk_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/fulldisk_60_temp.list','/whirlwind/goes17/wvc/goes17_fulldisk_loop_60.list')


# Fig for b/w fulldisk
fig2 = plt.figure(figsize=(20.0,20.0),dpi=80.)    
ax2 = fig2.add_subplot(1, 1, 1, projection=proj)
ax2.outline_patch.set_edgecolor('none')
im2 = ax2.imshow(a[:], extent=(xa[0],xa[-1],ya[-1],ya[0]), origin='upper', cmap=my_cmap2, vmin=162., vmax=330.)
ax2.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax2.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax2.add_feature(state_boundaries, linestyle=':')
ax2.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for fulldisk
#cbaxes2 = fig2.add_axes([0.135,0.12,0.755,0.02])
cbaxes2 = fig2.add_axes([0.135,0.12,0.255,0.005])
cbar2 = fig2.colorbar(im2, cax=cbaxes2, orientation='horizontal')
font_size = 10
#cbar2.set_label('Brightness Temperature (K)',size=18)
cbar2.ax.tick_params(labelsize=font_size, labelcolor='black')
cbar2.ax.xaxis.set_ticks_position('top')
cbar2.ax.xaxis.set_label_position('top')
text2 = ax2.text(0.005, 0.92, time_string,
    horizontalalignment='left', transform = ax2.transAxes,
    color='yellow', fontsize='14', weight='bold')

text2.set_path_effects(outline_effect)

filename2="/whirlwind/goes17/wv/fulldisk/"+ctime_file_string+"_fulldisk.jpg"
#filename2="/dev/shm/wv_"+ctime_file_string+"_fulldisk.jpg"
fig2.figimage(aoslogo,  0, fig2.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig2.savefig(filename2, bbox_inches='tight', pad_inches=0)
shutil.copy(filename2, "/whirlwind/goes17/wv/fulldisk/latest_fulldisk_1.jpg")
plt.close(fig2)

# FULLDISK B/W
file_list = glob.glob('/whirlwind/goes17/wv/fulldisk/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wv/fulldisk_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/fulldisk_24_temp.list','/whirlwind/goes17/wv/goes17_fulldisk_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wv/fulldisk_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/fulldisk_60_temp.list','/whirlwind/goes17/wv/goes17_fulldisk_loop_60.list')


# END Fulldisk

# START NPAC
# Fig for color NPAC region
fig3 = plt.figure(figsize=(npac_image_size_x,npac_image_size_y),dpi=80.)
ax3 = fig3.add_subplot(1, 1, 1, projection=proj)
ax3.outline_patch.set_edgecolor('none')
im3 = ax3.imshow(a[npac_image_crop_top:npac_image_crop_bottom,npac_image_crop_left:npac_image_crop_right], extent=(xa[npac_image_crop_left],xa[npac_image_crop_right],ya[npac_image_crop_bottom],ya[npac_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
ax3.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax3.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax3.add_feature(state_boundaries, linestyle=':')
ax3.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for npac
cbaxes3 = fig3.add_axes([0.135,0.14,0.755,0.02])
cbar3 = fig3.colorbar(im3, cax=cbaxes3, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar3.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar3.ax.xaxis.set_ticks_position('top')
cbar3.ax.xaxis.set_label_position('top')
text3 = ax3.text(0.005, 0.89, time_string,
    horizontalalignment='left', transform = ax3.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text3.set_path_effects(outline_effect)
#
filename3="/whirlwind/goes17/wvc/npac/"+ctime_file_string+"_npac.jpg"
#filename3="/dev/shm/wvc_"+ctime_file_string+"_npac.jpg"
fig3.figimage(aoslogo,  0, fig3.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig3.savefig(filename3, bbox_inches='tight', pad_inches=0)
shutil.copy(filename3, "/whirlwind/goes17/wvc/npac/latest_npac_1.jpg")
plt.close(fig3)

# NPAC Color
file_list = glob.glob('/whirlwind/goes17/wvc/npac/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wvc/npac_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/npac_24_temp.list','/whirlwind/goes17/wvc/goes17_npac_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wvc/npac_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/npac_60_temp.list','/whirlwind/goes17/wvc/goes17_npac_loop_60.list')


# Fig for b/w NPAC region
fig3 = plt.figure(figsize=(npac_image_size_x,npac_image_size_y),dpi=80.)
ax3 = fig3.add_subplot(1, 1, 1, projection=proj)
ax3.outline_patch.set_edgecolor('none')
im3 = ax3.imshow(a[npac_image_crop_top:npac_image_crop_bottom,npac_image_crop_left:npac_image_crop_right], extent=(xa[npac_image_crop_left],xa[npac_image_crop_right],ya[npac_image_crop_bottom],ya[npac_image_crop_top]), origin='upper',cmap=my_cmap2, vmin=162., vmax=330.)
ax3.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax3.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax3.add_feature(state_boundaries, linestyle=':')
ax3.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for npac
cbaxes3 = fig3.add_axes([0.135,0.14,0.755,0.02])
cbar3 = fig3.colorbar(im3, cax=cbaxes3, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar3.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar3.ax.xaxis.set_ticks_position('top')
cbar3.ax.xaxis.set_label_position('top')
text3 = ax3.text(0.005, 0.89, time_string,
    horizontalalignment='left', transform = ax3.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text3.set_path_effects(outline_effect)
#
filename3="/whirlwind/goes17/wv/npac/"+ctime_file_string+"_npac.jpg"
#filename3="/dev/shm/wv_"+ctime_file_string+"_npac.jpg"
fig3.figimage(aoslogo,  0, fig3.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig3.savefig(filename3, bbox_inches='tight', pad_inches=0)
shutil.copy(filename3, "/whirlwind/goes17/wv/npac/latest_npac_1.jpg")
plt.close(fig3)

# NPAC B/W
file_list = glob.glob('/whirlwind/goes17/wv/npac/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wv/npac_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/npac_24_temp.list','/whirlwind/goes17/wv/goes17_npac_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wv/npac_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/npac_60_temp.list','/whirlwind/goes17/wv/goes17_npac_loop_60.list')


# END NPAC

# START AK
# Fig for color AK region
fig4 = plt.figure(figsize=(ak_image_size_x,ak_image_size_y),dpi=80.)
ax4 = fig4.add_subplot(1, 1, 1, projection=proj)
ax4.outline_patch.set_edgecolor('none')
im4 = ax4.imshow(a[ak_image_crop_top:ak_image_crop_bottom,ak_image_crop_left:ak_image_crop_right], extent=(xa[ak_image_crop_left],xa[ak_image_crop_right],ya[ak_image_crop_bottom],ya[ak_image_crop_top]), origin='upper',cmap=my_cmap, vmin=162., vmax=330.)
ax4.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax4.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax4.add_feature(state_boundaries, linestyle=':')
ax4.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for ak
cbaxes4 = fig4.add_axes([0.135,0.14,0.755,0.02])
cbar4 = fig4.colorbar(im4, cax=cbaxes4, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar4.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar4.ax.xaxis.set_ticks_position('top')
cbar4.ax.xaxis.set_label_position('top')
text4 = ax4.text(0.005, 0.80, time_string,
    horizontalalignment='left', transform = ax4.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text4.set_path_effects(outline_effect)
#
filename4="/whirlwind/goes17/wvc/ak/"+ctime_file_string+"_ak.jpg"
#filename4="/dev/shm/wvc_"+ctime_file_string+"_ak.jpg"
fig4.figimage(aoslogo,  0, fig4.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig4.savefig(filename4, bbox_inches='tight', pad_inches=0)
shutil.copy(filename4, "/whirlwind/goes17/wvc/ak/latest_ak_1.jpg")

plt.close(fig4)

# AK Color
file_list = glob.glob('/whirlwind/goes17/wvc/ak/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wvc/ak_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/ak_24_temp.list','/whirlwind/goes17/wvc/goes17_ak_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wvc/ak_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wvc/ak_60_temp.list','/whirlwind/goes17/wvc/goes17_ak_loop_60.list')


# Fig for b/w AK region
fig4 = plt.figure(figsize=(ak_image_size_x,ak_image_size_y),dpi=80.)
ax4 = fig4.add_subplot(1, 1, 1, projection=proj)
ax4.outline_patch.set_edgecolor('none')
im4 = ax4.imshow(a[ak_image_crop_top:ak_image_crop_bottom,ak_image_crop_left:ak_image_crop_right], extent=(xa[ak_image_crop_left],xa[ak_image_crop_right],ya[ak_image_crop_bottom],ya[ak_image_crop_top]), origin='upper',cmap=my_cmap2, vmin=162., vmax=330.)
ax4.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax4.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax4.add_feature(state_boundaries, linestyle=':')
ax4.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for ak
cbaxes4 = fig4.add_axes([0.135,0.14,0.755,0.02])
cbar4 = fig4.colorbar(im4, cax=cbaxes4, orientation='horizontal')
font_size = 14
#cbar1.set_label('Brightness Temperature (K)',size=18)
cbar4.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar4.ax.xaxis.set_ticks_position('top')
cbar4.ax.xaxis.set_label_position('top')
text4 = ax4.text(0.005, 0.80, time_string,
    horizontalalignment='left', transform = ax4.transAxes,
    color='yellow', fontsize='14', weight='bold')
#
text4.set_path_effects(outline_effect)
#
filename4="/whirlwind/goes17/wv/ak/"+ctime_file_string+"_ak.jpg"
#filename4="/dev/shm/wv_"+ctime_file_string+"_ak.jpg"
fig4.figimage(aoslogo,  0, fig4.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig4.savefig(filename4, bbox_inches='tight', pad_inches=0)
shutil.copy(filename4, "/whirlwind/goes17/wv/ak/latest_ak_1.jpg")

plt.close(fig4)

# AK B/W
file_list = glob.glob('/whirlwind/goes17/wv/ak/2*.jpg')
file_list.sort()
#print("file list is ",file_list)

# 6h // 24 (36) images
thefile = open('/whirlwind/goes17/wv/ak_24_temp.list', 'w')
thelist = file_list[-36:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/ak_24_temp.list','/whirlwind/goes17/wv/goes17_ak_loop.list')

# 12h // 60 (72)images
thefile = open('/whirlwind/goes17/wv/ak_60_temp.list', 'w')
thelist = file_list[-72:]
#print ("thelist is ",thelist)

for item in thelist:
    head, tail = os.path.split(item)
    head, mid = os.path.split(head)
    thefile.write(mid + '/' + tail + '\n')
thefile.close
os.rename('/whirlwind/goes17/wv/ak_60_temp.list','/whirlwind/goes17/wv/goes17_ak_loop_60.list')


# END AK

# START full res full
# Fig for color full res fulldisk
fig9 = plt.figure(figsize=(image_columns/78.,image_rows/78.))
ax9 = fig9.add_subplot(1, 1, 1, projection=proj)
ax9.outline_patch.set_edgecolor('none')
im9 = ax9.imshow(a[:], extent=(xa[0],xa[-1],ya[-1],ya[0]), origin='upper', cmap=my_cmap, vmin=162., vmax=330.)
ax9.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax9.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax9.add_feature(state_boundaries, linestyle=':')
ax9.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for full
cbaxes9 = fig9.add_axes([0.135,0.15,0.755,0.02])
cbar9 = fig9.colorbar(im9, cax=cbaxes9, orientation='horizontal')
font_size = 18
#cbar9.set_label('Brightness Temperature (K)',size=18)
cbar9.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar9.ax.xaxis.set_ticks_position('top')
cbar9.ax.xaxis.set_label_position('top')
text9 = ax9.text(0.50, 0.97, time_string,
    horizontalalignment='center', transform = ax9.transAxes,
    color='yellow', fontsize='large', weight='bold')

text9.set_path_effects(outline_effect)
filename9="/whirlwind/goes17/wvc/fulldisk_full/"+ctime_file_string+"_fulldisk_full.jpg"
#filename9="/dev/shm/wvc_"+ctime_file_string+"_fulldisk_full.jpg"
fig9.figimage(aoslogo,  0, fig9.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig9.savefig(filename9, bbox_inches='tight', pad_inches=0)
shutil.copy(filename9, "/whirlwind/goes17/wvc/fulldisk_full/latest_fulldisk_full_1.jpg")
plt.close(fig9)

# Fig for b/w full res fulldisk
fig9 = plt.figure(figsize=(image_columns/78.,image_rows/78.))
ax9 = fig9.add_subplot(1, 1, 1, projection=proj)
ax9.outline_patch.set_edgecolor('none')
im9 = ax9.imshow(a[:], extent=(xa[0],xa[-1],ya[-1],ya[0]), origin='upper', cmap=my_cmap2, vmin=162., vmax=330.)
ax9.coastlines(resolution='50m', color='green')
# Add country borders with a thick line.
ax9.add_feature(cfeat.BORDERS, linewidth=1, edgecolor='green')
# Add the feature with dotted lines, denoted by ':'
ax9.add_feature(state_boundaries, linestyle=':')
ax9.add_feature(lat_lon_lines_15deg, linestyle=':')
# axes for full
cbaxes9 = fig9.add_axes([0.135,0.15,0.755,0.02])
cbar9 = fig9.colorbar(im9, cax=cbaxes9, orientation='horizontal')
font_size = 18
#cbar9.set_label('Brightness Temperature (K)',size=18)
cbar9.ax.tick_params(labelsize=font_size, labelcolor='yellow')
cbar9.ax.xaxis.set_ticks_position('top')
cbar9.ax.xaxis.set_label_position('top')
text9 = ax9.text(0.50, 0.97, time_string,
    horizontalalignment='center', transform = ax9.transAxes,
    color='yellow', fontsize='large', weight='bold')

text9.set_path_effects(outline_effect)
filename9="/whirlwind/goes17/wv/fulldisk_full/"+ctime_file_string+"_fulldisk_full.jpg"
#filename9="/dev/shm/wv_"+ctime_file_string+"_fulldisk_full.jpg"
fig9.figimage(aoslogo,  0, fig9.bbox.ymax*0.77 - aoslogoheight  , zorder=10)
fig9.savefig(filename9, bbox_inches='tight', pad_inches=0)
shutil.copy(filename9, "/whirlwind/goes17/wv/fulldisk_full/latest_fulldisk_full_1.jpg")
plt.close(fig9)

# END full res full

quit()

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
    
silentremove("/whirlwind/goes17/wvc/namer/latest_namer_24.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_23.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_24.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_22.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_23.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_21.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_22.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_20.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_21.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_19.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_20.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_18.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_19.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_17.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_18.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_16.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_17.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_15.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_16.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_14.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_15.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_13.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_14.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_12.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_13.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_11.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_12.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_10.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_11.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_9.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_10.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_8.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_9.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_7.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_8.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_6.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_7.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_5.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_6.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_4.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_5.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_3.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_4.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_2.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_3.jpg")
silentrename("/whirlwind/goes17/wvc/namer/latest_namer_1.jpg", "/whirlwind/goes17/wvc/namer/latest_namer_2.jpg")
    
shutil.copy(filename1, "/whirlwind/goes17/wvc/namer/latest_namer_1.jpg")


silentremove("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_24.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_23.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_24.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_22.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_23.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_21.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_22.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_20.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_21.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_19.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_20.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_18.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_19.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_17.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_18.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_16.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_17.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_15.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_16.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_14.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_15.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_13.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_14.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_12.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_13.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_11.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_12.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_10.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_11.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_9.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_10.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_8.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_9.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_7.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_8.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_6.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_7.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_5.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_6.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_4.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_5.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_3.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_4.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_2.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_3.jpg")
silentrename("/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_1.jpg", "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_2.jpg")
    
shutil.copy(filename2, "/whirlwind/goes17/wvc/fulldisk/latest_fulldisk_1.jpg")

silentremove("/whirlwind/goes17/wvc/npac/latest_npac_24.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_23.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_24.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_22.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_23.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_21.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_22.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_20.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_21.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_19.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_20.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_18.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_19.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_17.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_18.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_16.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_17.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_15.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_16.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_14.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_15.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_13.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_14.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_12.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_13.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_11.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_12.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_10.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_11.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_9.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_10.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_8.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_9.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_7.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_8.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_6.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_7.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_5.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_6.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_4.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_5.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_3.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_4.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_2.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_3.jpg")
silentrename("/whirlwind/goes17/wvc/npac/latest_npac_1.jpg", "/whirlwind/goes17/wvc/npac/latest_npac_2.jpg")
    
shutil.copy(filename3, "/whirlwind/goes17/wvc/npac/latest_npac_1.jpg")

silentremove("/whirlwind/goes17/wvc/ak/latest_ak_24.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_23.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_24.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_22.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_23.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_21.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_22.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_20.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_21.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_19.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_20.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_18.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_19.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_17.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_18.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_16.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_17.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_15.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_16.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_14.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_15.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_13.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_14.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_12.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_13.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_11.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_12.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_10.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_11.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_9.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_10.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_8.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_9.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_7.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_8.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_6.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_7.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_5.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_6.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_4.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_5.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_3.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_4.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_2.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_3.jpg")
silentrename("/whirlwind/goes17/wvc/ak/latest_ak_1.jpg", "/whirlwind/goes17/wvc/ak/latest_ak_2.jpg")
    
shutil.copy(filename4, "/whirlwind/goes17/wvc/ak/latest_ak_1.jpg")

shutil.copy(filename9, "/whirlwind/goes17/wvc/fulldisk_full/latest_fulldisk_full_1.jpg")

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


import sys
import pygrib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.basemap import Basemap
from scipy import interpolate
from scipy.interpolate import griddata
import numpy as np
#from pykrige.ok import OrdinaryKriging

def plotLayout(m):
    m.drawcoastlines()
    m.drawmapboundary()
    m.drawparallels(np.arange(-90.,120.,30.),labels=[1,0,0,0])
    m.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1])

if len(sys.argv) < 3:
    print 'Error: Missing arguments to script.'
    print 'Error: GRIB file name and the desired message record number should be given.'
    quit()

print 'Opening GRIB file ' + sys.argv[1] + ':'
grbs=pygrib.open(sys.argv[1])

msg = grbs[int(sys.argv[2])]
print 'Message corresponding to record number ' + sys.argv[2] + ':'
print msg

data=msg.values
print 'Variable values dimension (rows/cols):'
print data.shape
print 'Variable values Max/Min:'
print np.amax(data),np.amin(data)

lat,lon = msg.latlons() # Set the names of the latitude and longitude variables in the input GRIB file
print 'Latitude Max/Min:'
print np.amax(lat),np.amin(lat)
print 'Longitude Max/Min:'
print np.amax(lon),np.amin(lon)

print 'Grid resolution:'
resolutions = (abs(lat[0,0] - lat[1,0]), abs(lon[0,0] - lon[0,1]))
print resolutions

# We will use a coarser 5 x 5 degrees grid on which we do the interpolations
coarseDegree = 5.0
lat_coarse, lon_coarse = np.mgrid[np.amax(lat):np.amin(lat):36j, np.amin(lon):np.amax(lon):72j]

m=Basemap(projection='cyl', llcrnrlon=lon.min(), \
  urcrnrlon=lon.max(),llcrnrlat=lat.min(),urcrnrlat=lat.max(), \
  resolution='c')

# Original grid 
x, y = m(lon,lat)

# Coarser grid
x_coarse, y_coarse = m(lon_coarse, lat_coarse)

# Re-sample the data to the coarser grid
xrows_c = x_coarse.shape[0]
xcols_c = x_coarse.shape[1]
xy_coarse = []
data_coarse_list = []
for i in range(0, xrows_c):
    for j in range(0, xcols_c):
       xy_coarse.append((x_coarse[i,j],y_coarse[i,j]))

xrows = x.shape[0]
xcols = x.shape[1]

step1 = 0.0
step2 = 0.0
data_coarse_list.append(data[0,0])
for j in range(1, xcols):
    step2 = step2 + resolutions[1]
    if step2 == coarseDegree:
        data_coarse_list.append(data[0,j]) 
        step2 = 0.0

for i in range(1, xrows-1):
    step1 = step1 + resolutions[1]
    if step1 == coarseDegree:
       step2 = 0.0
       data_coarse_list.append(data[i,0])
       for j in range(1, xcols):
         step2 = step2 + resolutions[1]
         if step2 == coarseDegree:
            data_coarse_list.append(data[i,j])
            step2 = 0.0
       step1 = 0.0

data_coarse = np.zeros(shape=(xrows_c,xcols_c))
k = 0
for i in range(0, xrows_c):
    for j in range(0, xcols_c):
         data_coarse[i,j] = data_coarse_list[k]
         k = k + 1

# ..... Interpolations
# ... Linear and cubic methods
data_inperpol_linear  = griddata(xy_coarse, data_coarse_list, (x, y), method='linear')
data_inperpol_cubic   = griddata(xy_coarse, data_coarse_list, (x, y), method='cubic')

# ... Two-dimensional spline
tck = interpolate.bisplrep(x_coarse, y_coarse, data_coarse)

x_list = []
y_list = []
for i in range(0, xcols):
    x_list.append(x[0,i])
for i in range(0, xrows):
    y_list.append(y[i,0])

y_list.reverse() # bisplev requires the data to be in increasing order
z = interpolate.bisplev(x_list, y_list, tck)
z_transposed = np.transpose(z)
data_inperpol_spline = z_transposed[::-1] # Reverse the array

# ... Ordinary kriging method
# (Kriging was considered with the PyKrige toolkit but it was not working reliably with large grid size)
# OK = OrdinaryKriging(x_coarse_list, y_coarse_list, data_coarse_list, variogram_model='gaussian',
#                     verbose=False, enable_plotting=False)
# z = OK.execute('grid', x, y, backend='C')

# Plot original data
plt.figure(1,figsize=(15,11))
cs = m.pcolormesh(x,y,data,shading='flat',cmap=plt.cm.hot) #, norm=colors.LogNorm())

plotLayout(m)
plt.colorbar(cs,orientation='vertical', shrink=0.5)
plt.title('Temperature (K), isobaric level 100 hPa. Original grid.') 
#plt.savefig('TempOriginal.png') # Set the output file name
plt.show(block=False) # Show the plot in non-blocking mode 

# Plot the re-sampled coarser data
plt.figure(2,figsize=(15,11))
cs = m.pcolormesh(x_coarse,y_coarse,data_coarse,shading='flat',cmap=plt.cm.hot) #, norm=colors.LogNorm())

plotLayout(m)
m.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1])
m.drawmeridians(np.arange(0.,360.,60.),labels=[0,0,0,1])
plt.colorbar(cs,orientation='vertical', shrink=0.5)
plt.title('Temperature (K), isobaric level 100 hPa. Coarse grid.')
plt.show(block=False) # Show the plot in non-blocking mode 

# Plot interpolated data to original 0.25 x 0.25 degrees resolution grid with linear method
plt.figure(3,figsize=(15,11))
cs = m.pcolormesh(x,y,data_inperpol_linear,shading='flat',cmap=plt.cm.hot)

plotLayout(m)
plt.colorbar(cs,orientation='vertical', shrink=0.5)
plt.title('Temperature (K), isobaric level 100 hPa. Linear interpolation.')
plt.show(block=False)

# Plot interpolated data with cubic method
plt.figure(4,figsize=(15,11))
cs = m.pcolormesh(x,y,data_inperpol_cubic,shading='flat',cmap=plt.cm.hot)

plotLayout(m)
plt.colorbar(cs,orientation='vertical', shrink=0.5)
plt.title('Temperature (K), isobaric level 100 hPa. Cubic interpolation.')
plt.show(block=False)

# Plot interpolated data with spline method
plt.figure(5,figsize=(15,11))
cs = m.pcolormesh(x,y,data_inperpol_spline,shading='flat',cmap=plt.cm.hot)

plotLayout(m)
plt.colorbar(cs,orientation='vertical', shrink=0.5)
plt.title('Temperature (K), isobaric level 100 hPa. Spline interpolation.')
plt.show(block=False)

# Block main thread until all plots are closed
plt.show()


import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from matplotlib.backends.backend_pdf import PdfPages
import os
import csv
import re
import os
from astropy.io.ascii import read
from astropy.wcs import WCS
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
from matplotlib.patches import Rectangle, Circle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits

import funstools
from funstools import Cube2map, get_rms


from astropy import constants as co
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.visualization import simple_norm
from astropy.visualization.stretch import SinhStretch, LinearStretch, SqrtStretch, HistEqStretch
from astropy.visualization import ImageNormalize, MinMaxInterval

fdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'
funs_c18o = fdir+'SerB_A_C18O_match_cube.fits'
cube = Cube2map(funs_c18o, getrms='both', rmssize=300, velocity_smo=1, spatial_smo=1)
m0 = cube.moment0(vr=[-5,21])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(m0, origin='lower')


import pickle
from gausspyplus.plotting import pickle_load_file
pickle_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/sn2.5d2s1/output/gpy_training/SerB_A_C18O-training_set_300_fwhm50_s2.5d2s1.pickle'
data = pickle_load_file(pickle_file)
for i in range(len(data['location'])):
    y = data['location'][i][0]
    x = data['location'][i][1]
    ax.plot(x,y,color='green', marker='o', markersize=3)

plt.title('SerB Cluster A: Selected 300 pixels for the Training Set', fontsize=9)
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/sn2.5d2s1/output/gpy_training/SerB_A_C18O-training_set_300_fwhm50_s2.5d2.pdf')


####---------------------------

fdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
funs_c18o = fdir+'SerB_C18O_match_cube.fits'
cube = Cube2map(funs_c18o, getrms='both', rmssize=300, velocity_smo=1, spatial_smo=1)
m0 = cube.moment0(vr=[-5,21])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(m0, origin='lower')


import pickle
from gausspyplus.plotting import pickle_load_file
pickle_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o/c18o_m521_training/SerB_C18O_m521-training_set_500_spectra.pickle'
data = pickle_load_file(pickle_file)
for i in range(len(data['location'])):
    y = data['location'][i][0]
    x = data['location'][i][1]
    ax.plot(x,y,color='green', marker='o', markersize=3)

plt.title('SerB: Selected 500 pixels for the Training Set', fontsize=9)
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o/c18o_m521_training/SerB_C18O_m521-training_set_500_spectra_position.pdf')



####---------------------------

fdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
funs_c18o = fdir+'SerB_C18O_match_cube.fits'
cube = Cube2map(funs_c18o, getrms='both', rmssize=300, velocity_smo=1, spatial_smo=1)
m0 = cube.moment0(vr=[-5,21])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(m0, origin='lower')


import pickle
from gausspyplus.plotting import pickle_load_file
pickle_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o/c18o_training/SerB_C18O-training_set_150_spectra_ord2_mfwhm35.pickle'
data = pickle_load_file(pickle_file)
for i in range(len(data['location'])):
    y = data['location'][i][0]
    x = data['location'][i][1]
    ax.plot(x,y,color='green', marker='o', markersize=3)

plt.title('SerB: Selected 150 pixels for the Training Set', fontsize=9)
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o/c18o_training/SerB_C18O-training_set_150_spectra_ord2_mfwhm35_position.pdf')


####---------------------------------

fdir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/CropCubeAB/FUNS/'
funs_c18o = fdir+'SerB_A_C18O_match_cube.fits'
cube = Cube2map(funs_c18o, getrms='both', rmssize=300, velocity_smo=1, spatial_smo=1)
m0 = cube.moment0(vr=[-5,21])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(m0, origin='lower')


import pickle
from gausspyplus.plotting import pickle_load_file
pickle_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_training/SerB_A_C18O-training_set_sn3n3_fwhm25_300_spectra.pickle'
data = pickle_load_file(pickle_file)
for i in range(len(data['location'])):
    y = data['location'][i][0]
    x = data['location'][i][1]
    ax.plot(x,y,color='green', marker='o', markersize=3)

plt.title('SerB Cluster A: Selected 300 pixels for the Training Set (sn3n3 fwhm25)', fontsize=9)
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_training/SerB_A_C18O-training_set_sn3n3_fwhm25_300_spectra_position.pdf')




####-----------
import pickle
from gausspyplus.plotting import pickle_load_file
pickle_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o/gpy_prepared/SerB_C18O_match_cube.pickle'
data = pickle_load_file(pickle_file)




#------------- after decomposition
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

#----sn3n3 fwhm25
filepath = os.path.join('decomposition_c18o_A_sn3n3', 'gpy_maps', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_rchi2_map.fits')
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 12
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('decomposition_c18o_A_sn3n3', 'gpy_maps', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_component_map.fits'))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After first decomposition with improved fitting routine (Cluster A sn3n3 fwhm25)')

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_plots/SerB_A_C18O_sn3n3_fwhm25_rch2_Nc_g+_fit_fin.pdf')



#----sn3n3 no max_fwhm constraint
filepath = os.path.join('decomposition_c18o_A_sn3n3', 'gpy_maps', 'sn3n3', 'SerB_A_C18O_match_cube_g+_rchi2_map.fits')
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 12
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('decomposition_c18o_A_sn3n3', 'gpy_maps', 'sn3n3','SerB_A_C18O_match_cube_g+_component_map.fits'))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After first decomposition with improved fitting routine (Cluster A sn3n3)')

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_plots/SerB_A_C18O_sn3n3_rch2_Nc_g+_fit_fin.pdf')







#------------- after phase 1
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

#--fwhm25
filepath = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_fit_fin_sf-p1_rchi2_map.fits')
rchi2 = fits.getdata(filepath)
header = fits.getheader(filepath)
wcs = WCS(header)

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2.5
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_fit_fin_sf-p1_component_map.fits'))

vmax = 6
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 1 (Cluster A sn3n3 fwhm25)')

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_plots/SerB_A_C18O_sn3n3_fwhm25_rch2_Nc_g+_fit_fin_sf-p1.pdf')



#--sn3n3

filepath = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps/sn3n3', 'SerB_A_C18O_match_cube_sn3n3_g+_fit_fin_sf-p1_rchi2_map.fits')
rchi2 = fits.getdata(filepath)
header = fits.getheader(filepath)
wcs = WCS(header)

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2.5
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps/sn3n3', 'SerB_A_C18O_match_cube_sn3n3_g+_fit_fin_sf-p1_component_map.fits'))

vmax = 6
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 1 (Cluster A sn3n3)')

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_plots/SerB_A_C18O_sn3n3_rch2_Nc_g+_fit_fin_sf-p1.pdf')




#------------- after phase 2
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

#--fwhm25
filepath = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_fit_fin_sf-p2_rchi2_map.fits')
rchi2 = fits.getdata(filepath)
header = fits.getheader(filepath)
wcs = WCS(header)

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2.5
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_fit_fin_sf-p2_component_map.fits'))

vmax = 6
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 2 (Cluster A sn3n3 fwhm25)')

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_plots/SerB_A_C18O_sn3n3_fwhm25_rch2_Nc_g+_fit_fin_sf-p2.pdf')



#--sn3n3

filepath = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps/sn3n3', 'SerB_A_C18O_match_cube_sn3n3_g+_fit_fin_sf-p2_rchi2_map.fits')
rchi2 = fits.getdata(filepath)
header = fits.getheader(filepath)
wcs = WCS(header)

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2.5
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_maps/sn3n3', 'SerB_A_C18O_match_cube_sn3n3_g+_fit_fin_sf-p2_component_map.fits'))

vmax = 6
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 2 (Cluster A sn3n3)')

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A_sn3n3/gpy_plots/SerB_A_C18O_sn3n3_rch2_Nc_g+_fit_fin_sf-p2.pdf')







###### 7/26 8 sets test results comparision


#------------- after decomposition
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

#----decomposed
versionname = 'sn3.0d3s2'

filepath = os.path.join('output', 'gpy_maps', 'SerB_A_C18O_match_cube_{}_g+_rchi2_map.fits'.format(versionname))
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('output', 'gpy_maps', 'SerB_A_C18O_match_cube_{}_g+_component_map.fits'.format(versionname)))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After first decomposition with improved fitting routine (Cluster A {})'.format(versionname))

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/{}/output/gpy_plots/SerB_A_C18O_match_cube_{}_rchi2_Nc_g+_fit_fin.pdf'.format(versionname,versionname))



#------------- after phase 1
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')


#---- after refitting phase 1
versionname = 'sn3.0d3s2'

filepath = os.path.join('output', 'gpy_maps', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p1_rchi2_map.fits'.format(versionname))
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('output', 'gpy_maps', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p1_component_map.fits'.format(versionname)))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 1 (Cluster A {})'.format(versionname))

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/{}/output/gpy_plots/SerB_A_C18O_match_cube_{}_rchi2_Nc_g+_fit_fin_sf-p1.pdf'.format(versionname,versionname))




#------------- after phase 2
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')


#---- after refitting phase 2
versionname = 'sn3.0d3s2'

filepath = os.path.join('output', 'gpy_maps', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_rchi2_map.fits'.format(versionname))
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('output', 'gpy_maps', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2_component_map.fits'.format(versionname)))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 2 (Cluster A {})'.format(versionname))

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/{}/output/gpy_plots/SerB_A_C18O_match_cube_{}_rchi2_Nc_g+_fit_fin_sf-p2.pdf'.format(versionname, versionname))





###### 8/07 SNR=5, ord=3, significance=5


#------------- after decomposition
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')

#----decomposed
versionname = 'sn5.0d3s1'

filepath = os.path.join('output', 'gpy_maps', 'SerB_C18O_match_cube_{}_g+_rchi2_map.fits'.format(versionname))
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('output', 'gpy_maps', 'SerB_C18O_match_cube_{}_g+_component_map.fits'.format(versionname)))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After first decomposition with improved fitting routine (Main {})'.format(versionname))

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/{}/output/gpy_plots/SerB_C18O_match_cube_{}_rchi2_Nc_g+_fit_fin.pdf'.format(versionname,versionname))



#------------- after phase 1
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')


#---- after refitting phase 1
versionname = 'sn5.0d3s1'

filepath = os.path.join('output', 'gpy_maps', 'SerB_C18O_match_cube_{}_g+_fit_fin_sf-p1_rchi2_map.fits'.format(versionname))
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('output', 'gpy_maps', 'SerB_C18O_match_cube_{}_g+_fit_fin_sf-p1_component_map.fits'.format(versionname)))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 1 (SerB {})'.format(versionname))

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/{}/output/gpy_plots/SerB_C18O_match_cube_{}_rchi2_Nc_g+_fit_fin_sf-p1.pdf'.format(versionname,versionname))




#------------- after phase 2
import os

from astropy.io import fits

import matplotlib
import matplotlib.pyplot as plt
from pylab import cm

from astropy.wcs import WCS

from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap



def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)


def add_style(ax):
    ax.set_xlabel('RA')
    ax.set_ylabel('Dec')


#---- after refitting phase 2
versionname = 'sn5.0d3s1'

filepath = os.path.join('output', 'gpy_maps', 'SerB_C18O_match_cube_{}_g+_fit_fin_sf-p2_rchi2_map.fits'.format(versionname))
rchi2 = fits.getdata(filepath)
wcs = WCS(fits.getheader(filepath))

fig, axes = plt.subplots(ncols=2, figsize=(12, 4), subplot_kw=dict(projection=wcs))

ax = axes.flatten()[0]

vmin = min(rchi2.flatten())
vmax = 2
new_cmap = get_cmap_rchi2(vmin, vmax)

img_rchi2 = ax.imshow(rchi2, cmap=new_cmap, vmin=vmin, vmax=vmax)
fig.colorbar(img_rchi2, ax=ax, extend='max')
ax.set_title('$\chi_{\mathrm{red}}^{2}$ map')
add_style(ax)

ax = axes.flatten()[1]

ncomps = fits.getdata(os.path.join('output', 'gpy_maps', 'SerB_C18O_match_cube_{}_g+_fit_fin_sf-p2_component_map.fits'.format(versionname)))

vmax = 7
new_cmap = cm.get_cmap('Spectral_r', vmax + 1)

img_ncomps = ax.imshow(ncomps, cmap=new_cmap, vmin=0, vmax=vmax)
fig.colorbar(img_ncomps, ax=ax)
ax.set_title('Number of fitted components')
add_style(ax)

fig.suptitle('After spatially coherent refitting - phase 2 (Cluster A {})'.format(versionname))

plt.show()
fig.savefig('/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/{}/output/gpy_plots/SerB_C18O_match_cube_{}_rchi2_Nc_g+_fit_fin_sf-p2.pdf'.format(versionname, versionname))

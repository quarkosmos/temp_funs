"' to save a cropped data in velocity range to a fits file '"

import numpy as np
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
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits

import funstools
from funstools import Cube2map
#from run_funstools_Serp import get_filename
#from run_funstools_Serp import get_maxmin_maps

from astropy import constants as co
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.nddata import Cutout2D
from astropy.visualization import simple_norm

from spectral_cube import SpectralCube
import pyspeckit


path='/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/match_cube/'
cube_C18O_v10 = path+'SerB_C18O_v10_match_cube.fits'
cube_C18O_v06 = path+'SerB_C18O_match_cube.fits'

outpath = '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/vrcut_cube/'
cube_C18O_v10_m521 = outpath+'SerB_C18O_v10_m521_match_cube.fits'
cube_C18O_v06_m521 = outpath+'SerB_C18O_m521_match_cube.fits'
cube_C18O_v10_0313 = outpath+'SerB_C18O_v10_0313_match_cube.fits'
cube_C18O_v06_0313 = outpath+'SerB_C18O_0313_match_cube.fits'

funslist = [cube_C18O_v10, cube_C18O_v06]
v_start = [[-5, 3],[-5.02,2.96]]
wlist = [[cube_C18O_v10_m521,cube_C18O_v10_0313], [cube_C18O_v06_m521, cube_C18O_v06_0313]]

for i in range(len(funslist)):
    for j in range(len(v_start)):
        cube = Cube2map(funslist[i])
        #restfreq = 0.9317376370000E+11 * u.Hz # N2H+ 23 - 12 frequency
        #freq_to_vel = u.doppler_radio(restfreq)
        #freq_0112 = 93176.2522 * u.MHz # N2H+ isolated component
        #freq_0112_vel = (freq_0112).to(u.km / u.s, equivalencies=freq_to_vel)
        #flt_freq_0112_vel = float(freq_0112_vel / (u.km / u.s))
        lsr_vel_ch = cube.v2ch(8)
        #new_lsr_vel = 0.8000000000000E+04 + float(freq_0112_vel.to(u.m / u.s) / (u.m/ u.s)) # in m/s
        delch = lsr_vel_ch - cube.v2ch(v_start[i][j])
        if i == 0:
            data = cube.data[lsr_vel_ch-delch:lsr_vel_ch+delch, : , :]
        else:
            data = cube.data[lsr_vel_ch-delch:lsr_vel_ch+delch, : , :]
        h = cube.header
        if i == 0:
            h['ALTRPIX'] = delch+1 #0.5000000000000E+03  / Alternate freq/vel ref pixel
        else:
            h['ALTRPIX'] = delch+1
        #h['ALTRVAL'] = 0.931762522E+11  #/ Alternate freq/vel ref value
        #h['RESTFREQ'] = 0.931762522E+11 #/ Rest frequency (HZ)
        #h['VELO-LSR'] = new_lsr_vel #-0.00064E+04 #0.8000000000000E+04  / Velocity of reference channel (M/S)
        #h['VELREF'] = delch
        h['CRVAL3'] = v_start[i][j]
        h['COMMENT'] = 'Velocity Cropped map'
        fits.writeto(wlist[i][j], data=data, header=h, overwrite=True)






###################################################
"' temparary space to plot spectra from gpy+'"




path_to_file = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01', 'vrcut_cube', 'SerB_C18O_m521_match_cube.fits')
#  Directory to which all files produced by GaussPy+ will get saved.
dirpath_gpy = 'decomposition_c18o'
#  Number of spectra included in the training set. We recommend to have at least 250 spectra for a good training set.
filename_out = 'SerB_C18O_m521-training_set_{}_spectra.pickle'.format(500)

#  Filepath to pickled dictionary of the training set.
path_to_training_set = os.path.join(dirpath_gpy, 'c18o_m521_training', filename_out)
#  Directory in which the plots are saved.
path_to_plots = os.path.join(dirpath_gpy, 'c18o_m521_training')
plot_spectra(path_to_training_set,
             path_to_plots=path_to_plots,
             training_set=True,
             n_spectra=500  # Plot 20 random spectra of the training set.
             )


#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

dirpath_gpy = 'decomposition_c18o_A'
filename_out = 'SerB_A_C18O-training_set_{}_spectra.pickle'.format(250)
path_to_training_set = os.path.join(dirpath_gpy, 'c18o_A_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'c18o_A_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=450,
                 fig_max_chan=550, cols=3,
                 n_spectra=50  # Plot 20 random spectra of the training set.
                 )



#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

dirpath_gpy = 'decomposition_c18o_A'
filename_out = 'SerB_A_C18O-training_set_{}_spectra_ord3.pickle'.format(250)
path_to_training_set = os.path.join(dirpath_gpy, 'c18o_A_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'c18o_A_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=450,
                 fig_max_chan=550, cols=3,
                 n_spectra=50  # Plot 20 random spectra of the training set.
                 )


#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

dirpath_gpy = 'decomposition_c18o'
filename_out = 'SerB_C18O_m521-training_set_{}_spectra.pickle'.format(500) #ord 3
path_to_training_set = os.path.join(dirpath_gpy, 'c18o_m521_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'c18o_m521_training')
plot_spectra(path_to_training_set,
             path_to_plots=path_to_plots,
             training_set=True,fig_min_chan=118,
             fig_max_chan=318, cols=3,
             n_spectra=500  # Plot 20 random spectra of the training set.
             )


#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

dirpath_gpy = 'decomposition_c18o'
filename_out = 'SerB_C18O-training_set_{}_spectra_ord2_mfwhm35.pickle'.format(150)
path_to_training_set = os.path.join(dirpath_gpy, 'c18o_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'c18o_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=440,
                 fig_max_chan=560, cols=3,
                 n_spectra=150  # Plot 20 random spectra of the training set.
                 )



#------------
# to draw some pixels in pareparing set
import os
from gausspyplus.plotting import plot_spectra

#  Filepath to pickled dictionary of the prepared data.
path_to_pickled_file = os.path.join(
    'decomposition_c18o', 'gpy_prepared', 'SerB_C18O_match_cube.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [30, 39], 'y': [175, 184]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, cols=3)

#------------
# to draw some pixels in decomposed set
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o', 'gpy_prepared', 'SerB_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o', 'gpy_decomposed', 'SerB_C18O_match_cube_g+_fit_fin.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [55, 59], 'y': [60, 64]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)


#------------
# to draw some pixels in decomposed set
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o', 'gpy_prepared', 'SerB_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o', 'gpy_decomposed', 'SerB_C18O_match_cube_g+_fit_fin_sf-p1.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [55, 59], 'y': [60, 64]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)




#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A', 'gpy_decomposed', 'SerB_A_C18O_match_cube_g+_fit_fin.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [6, 10], 'y': [34, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)


#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A', 'gpy_decomposed', 'SerB_A_C18O_match_cube_g+_fit_fin_sf-p1.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [6, 10], 'y': [34, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)

#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A', 'gpy_decomposed', 'SerB_A_C18O_match_cube_g+_fit_fin_sf-p2.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [6, 10], 'y': [34, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)


#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A', 'gpy_decomposed', 'SerB_A_C18O_match_cube_g+_fit_fin_sf-p2_sf-p2.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [6, 10], 'y': [34, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)


#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A', 'gpy_decomposed', 'SerB_A_C18O_match_cube_g+_fit_fin_sf-p2_sf-p2_sf-p2.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [6, 10], 'y': [34, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)




#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A', 'gpy_decomposed', 'SerB_A_C18O_match_cube_g+_fit_fin_sf-p2_sf-p2_sf-p2_sf-p2.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [6, 10], 'y': [34, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)



#---------------------

#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_decomposed', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_fit_fin.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)



#------------
# to draw some pixels in decomposed set Cluster A fwhm25
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_decomposed', 'SerB_A_C18O_match_cube_sn3n3_fwhm25_g+_fit_fin_sf-p2.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [13, 15], 'y': [36, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)



#------------
# to draw some pixels in decomposed set Cluster A
import os
from gausspyplus.plotting import plot_spectra

path_to_pickled_file = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_prepared/sn3n3', 'SerB_A_C18O_match_cube_sn3n3.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_decomposed/sn3n3', 'SerB_A_C18O_match_cube_sn3n3_g+_fit_fin_sf-p2.pickle')
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'decomposition_c18o_A_sn3n3', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [13, 15], 'y': [36, 38]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560)



# sn2.5 training set order 비{}
#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

dirpath_gpy = 'output'
filename_out = 'SerB_C18O_A-training_set_{}_spectra_ord{}_fwhm{}.pickle'.format(30, 6, 50)
path_to_training_set = os.path.join(dirpath_gpy, 'gpy_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'gpy_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=440,
                 fig_max_chan=560, cols=3,
                 n_spectra=30  # Plot 20 random spectra of the training set.
                 )




# week of 7/26 8 sets:
# s2.5d2s1, s2.5d2s2, s2.5d3s1, s2.5d3s2
# s3.0d2s1, s3.0d2s2, s3.0d3s1, s3.0d3s2
#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

#step 1
dirpath_gpy = 'output'
filename_out = 'SerB_C18O_A-training_set_{}_fwhm{}_{}.pickle'.format(300, 50, 'sn3.0d2s1')
path_to_training_set = os.path.join(dirpath_gpy, 'gpy_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'gpy_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=440,
                 fig_max_chan=560, cols=3,
                 n_spectra=300, snr=3.0  # Plot 20 random spectra of the training set.
                 )

#step 2: train machine

#step 3: prepare data_a

#step 4: decompose
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin.pickle'.format('sn2.5d2s1'))
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)

#step 5: refit phase 1
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p1.pickle'.format('sn2.5d2s1'))
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)


#step 6: refit phase 2
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2.pickle'.format('sn2.5d2s1'))
#  Directory in which the plots are saved.n
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)





# week of 8/2 1 sets:
# s5.0d3s1
#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

#step 1
dirpath_gpy = 'output'
filename_out = 'SerB_C18O_A-training_set_{}_fwhm{}_{}.pickle'.format(300, 25, 'sn10.0d2s1')
path_to_training_set = os.path.join(dirpath_gpy, 'gpy_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'gpy_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=440,
                 fig_max_chan=560, cols=3,
                 n_spectra=166, snr=10.0  # Plot 20 random spectra of the training set.
                 )

#step 2: train machine

#step 3: prepare data_a

#step 4: decompose
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin.pickle'.format('sn5.0d3s1'))
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)

#step 5: refit phase 1
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p1.pickle'.format('sn5.0d3s1'))
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)


#step 6: refit phase 2
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2.pickle'.format('sn5.0d3s1'))
#  Directory in which the plots are saved.n
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)







# week of 8/7 1 sets:
# sn5.0d3s1 : significance=5
#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

#step 1
dirpath_gpy = 'output'
filename_out = 'SerB_C18O_A-training_set_{}_fwhm{}_{}.pickle'.format(500, 50, 'sn5.0d3s1')
path_to_training_set = os.path.join(dirpath_gpy, 'gpy_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'gpy_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=440,
                 fig_max_chan=560, cols=3,
                 n_spectra=500, snr=5.0  # Plot 20 random spectra of the training set.
                 )

#step 2: train machine

#step 3: prepare data_a

#step 4: decompose
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_C18O_match_cube_{}_g+_fit_fin.pickle'.format('sn5.0d3s1'))
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)

#step 5: refit phase 1
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p1.pickle'.format('sn5.0d3s1'))
#  Directory in which the plots are saved.
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)


#step 6: refit phase 2
path_to_pickled_file = os.path.join(
    'output', 'gpy_prepared', 'SerB_A_C18O_match_cube.pickle')
#  Filepath to pickled dictionary with the decomposition results
path_to_decomp_pickle = os.path.join(
    'output', 'gpy_decomposed', 'SerB_A_C18O_match_cube_{}_g+_fit_fin_sf-p2.pickle'.format('sn5.0d3s1'))
#  Directory in which the plots are saved.n
path_to_plots = os.path.join(
    'output', 'gpy_plots')
#  Here we select a subregion of the data cube, whose spectra we want to plot.
pixel_range = {'x': [10, 14], 'y': [30, 34]}
plot_spectra(path_to_pickled_file, path_to_plots=path_to_plots,
             path_to_decomp_pickle=path_to_decomp_pickle,
             signal_ranges=True, pixel_range=pixel_range,
             fig_min_chan=440,
             fig_max_chan=560, snr=2.5)





# week of 8/8 1 sets:
# sn3.0d2s1 : significance=5, max_fwhm=30
#------------------------
import os
from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra

#step 1
dirpath_gpy = 'output'
filename_out = 'SerB_C18O-training_set_{}_fwhm{}_{}.pickle'.format(500, 30, 'sn3.0d2s1')
path_to_training_set = os.path.join(dirpath_gpy, 'gpy_training', filename_out)
path_to_plots = os.path.join(dirpath_gpy, 'gpy_training')
plot_spectra(path_to_training_set,
                 path_to_plots=path_to_plots,
                 training_set=True,
                 fig_min_chan=440,
                 fig_max_chan=560, cols=3,
                 n_spectra=500, snr=3.0  # Plot 20 random spectra of the training set.
                 )


import pickle
from gausspyplus.plotting import pickle_load_file
pickle_file = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/decomposition_c18o_A/c18o_A_training/SerB_A_C18O-training_set_250_spectra.pickle'
data = pickle_load_file(pickle_file)
ymin, xmin = min(data['location'])

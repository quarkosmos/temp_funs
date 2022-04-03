# @Author: Kyounghee Kim
# @Date:   2021-07-06
# @Filename: training_set--c18o-A.py
# @Last modified for:   application to FUNS Serpens data (the original from GaussPy+)
# @Last modified time: 2021-07-23

import os

from gausspyplus.training_set import GaussPyTrainingSet
from gausspyplus.plotting import plot_spectra


def main():
    #  Initialize the 'GaussPyTrainingSet' class and read in the parameter settings from 'gausspy+.ini'.
    training = GaussPyTrainingSet(config_file='gausspy+sn5.0d2s1.ini')

    #  The following lines will override the corresponding parameter settings defined in 'gausspy+.ini'.

    #  Path to the FITS cube.
    training.path_to_file = os.path.join(
        '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01', 'rect_cube', 'SerB_C18O_match_v06_cube_rect.fits')
    #  Directory to which all files produced by GaussPy+ will get saved.
    training.dirpath_gpy = 'output'
    #  Number of spectra included in the training set. We recommend to have at least 250 spectra for a good training set.
    training.n_spectra = 500
    #  (Optional) The initial seed that is used to create pseudorandom numbers. Change this value in case the spectra chosen for the training set are not ideal.
    training.random_seed = 111
    #  (Optional) We set the upper limit for the reduced chi-square value to a lower number to only include good fits in the training sample
    training.rchi2_limit = 1.5
    #  (Optional) This will enforce a maximum upper limit for the FWHM value of fitted Gaussian components, in this case 50 channels. We recommended to use this upper limit for the FWHM only for the creation of the training set.
    training.max_fwhm = 67 #4.02km/s
    
    
    # by KHK
    #training.min_channels = 85 # ~5 km/s
    
    #  (Optional) Here we specify the filename for the resulting pickled dictionary file. If 'filename_out' is not supplied it will be automatically generated.
    training.filename_out = \
        'SerB_C18O_A-training_set_{}_fwhm{}_sn5.0d2s1.pickle'.format(training.n_spectra, training.max_fwhm)


    training.decompose_spectra()  # Create the training set.

if __name__ == "__main__":
    main()

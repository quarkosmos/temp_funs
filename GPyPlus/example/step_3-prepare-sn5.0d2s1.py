# @Author: riener
# @Date:   2019-04-02T17:23:06+02:00
# @Filename: prepare--grs.py
# @Last modified by:   riener
# @Last modified time: 2019-04-08T10:33:45+02:00


import os

from gausspyplus.prepare import GaussPyPrepare
from gausspyplus.plotting import plot_spectra


def main():
    #  Initialize the 'GaussPyPrepare' class and read in the parameter settings from 'gausspy+.ini'.
    prepare = GaussPyPrepare(config_file='gausspy+sn5.0d2s1.ini')

    #  The following lines will override the corresponding parameter settings defined in 'gausspy+.ini'.

    #  Path to the FITS cube.
    prepare.path_to_file = os.path.join(
        '/Users/khkim/Dropbox/Work5/FUNS-Serp/datared_01/rect_cube', 'SerB_C18O_match_v06_cube_rect.fits')
    #  Directory in which all files produced by GaussPy+ are saved.
    prepare.dirpath_gpy = 'output'
    #  Prepare the data cube for the decomposition
    prepare.prepare_cube()
    #  (Optional) Produce a FITS image with the estimated root-mean-square values
    prepare.produce_noise_map()



if __name__ == "__main__":
    main()



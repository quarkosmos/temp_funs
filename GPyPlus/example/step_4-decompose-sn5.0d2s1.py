# @Author: riener
# @Date:   2019-04-02T17:42:46+02:00
# @Filename: decompose--grs.py
# @Last modified by:   riener
# @Last modified time: 31-05-2019


import os

from gausspyplus.decompose import GaussPyDecompose
from gausspyplus.plotting import plot_spectra


def main():
    #  Initialize the 'GaussPyDecompose' class and read in the parameter settings from 'gausspy+.ini'.
    decompose = GaussPyDecompose(config_file='gausspy+sn5.0d2s1.ini')

    #  The following lines will override the corresponding parameter settings defined in 'gausspy+.ini'.

    #  Filepath to pickled dictionary of the prepared data.
    decompose.path_to_pickle_file = os.path.join(
        'output', 'gpy_prepared', 'SerB_C18O_vsm5.0_ssm3.0_v06_cube_rect.pickle')
    #  First smoothing parameter
    decompose.alpha1 = 3.12
    #  Second smoothing parameter
    decompose.alpha2 = 5.36
    #  Suffix for the filename of the pickled dictionary with the decomposition results.
    decompose.suffix = '_sn5.0d2s1_g+'
    
    
    #  Start the decomposition.
    decompose.decompose()

    #  (Optional) Produce a FITS image showing the number of fitted components
    decompose.produce_component_map()
    #  (Optional) Produce a FITS image showing the reduced chi-square values
    decompose.produce_rchi2_map()

if __name__ == "__main__":
    main()

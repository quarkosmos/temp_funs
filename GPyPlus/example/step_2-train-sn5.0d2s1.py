# @Author: riener
# @Date:   2019-04-02T16:40:08+02:00
# @Filename: train--grs.py
# @Last modified by:   riener
# @Last modified time: 2019-04-08T10:27:00+02:00


import os

from gausspyplus.training import GaussPyTraining


def main():
    #  Initialize the 'GaussPyTraining' class and read in the parameter settings from 'gausspy+.ini'.
    train = GaussPyTraining(config_file='gausspy+sn5.0d2s1.ini')

    #  The following lines will override the corresponding parameter settings defined in 'gausspy+.ini'.

    #  Directory in which all files produced by GaussPy+ are saved.
    train.dirpath_gpy = 'output'
    #train.gpy_dirpath_ = 'decomposition_c18o'
    #  Filepath to the training set.
    train.path_to_training_set = os.path.join(
        train.dirpath_gpy, 'gpy_training',
        'SerB_C18O_A-training_set_{}_fwhm67_sn5.0d2s1.pickle'.format(train.n_spectra))
    #  We select the two-phase-decomposition that uses two smoothing parameters.
    train.two_phase_decomposition = True
    #  Initial value for the first smoothing parameter.
    train.alpha1_initial = 2.
    #  Initial value for the second smoothing parameter.
    train.alpha2_initial = 6.
    #  Start the training.
    train.training()


if __name__ == "__main__":
    main()

# @Author: Manuel Riener <riener>
# @Date:   04-05-2020
# @Email:  riener@mpia-hd.mpg.de
# @Last modified by:   riener
# @Last modified time: 04-05-2020

import os

from gausspyplus.finalize import Finalize


def main():
    #  Initialize the 'Finalize' class and read in the parameter settings from 'gausspy+.ini'.
    finalize = Finalize(config_file='gausspy+sn5.0d2s1.ini')

    #  The following lines will override the corresponding parameter settings defined in 'gausspy+.ini'.

    #  filepath to the pickled dictionary of the prepared data
    finalize.path_to_pickle_file = os.path.join(
        'output', 'gpy_prepared', 'SerB_C18O_match_v06_cube_rect.pickle')
    #  Filepath to the pickled dictionary of the decomposition results
    finalize.path_to_decomp_file = os.path.join(
        'output', 'gpy_decomposed',
        'SerB_C18O_vsm5.0_ssm3.0_v06_cube_rect_sn5.0d2s1_g+_fit_fin_sf-p2.pickle')

    finalize.finalize_dct()
    finalize.make_table()


if __name__ == "__main__":
    main()

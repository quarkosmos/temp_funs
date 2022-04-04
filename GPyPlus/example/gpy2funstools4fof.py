
#############
## gpy finalized data  to the funstools decomposed data table format

from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn3.0d2s1' # the decomposed paramter option title (i.e., directory name)

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option,'output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_match_cube_{}_g+_fit_fin_sf-p1_finalized.dat'.format(dc_option)
gpy2fitres = 'SerB_C18O_match_cube_{}_g+_fit_fin_sf-p1_finalized_fitres.dat'.format(dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
#        print(j)
    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)

xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, snr, dv, int_tot], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'))
fr.write(os.path.join(path_to_decomp_data, gpy2fitres), format='ascii', overwrite='True')


############# smoothed data
## gpy finalized data  to the funstools decomposed data table format

from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn5.0d3s1' # the decomposed paramter option title (i.e., directory name)

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option+'_sv2521','output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_vsm5.0_ssm2.0_v06_cube_{}_g+_fit_fin_sf-p1_finalized.dat'.format(dc_option)
gpy2fitres = 'SerB_C18O_vsm5.0_ssm2.0_v06_cube_{}_g+_fit_fin_sf-p1_finalized_fitres.dat'.format(dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
#        print(j)

    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)

xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, snr, dv, int_tot], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'))
fr.write(os.path.join(path_to_decomp_data, gpy2fitres), format='ascii', overwrite='True')




############# smoothed data (step 4)
## gpy step 4 decomposed data  to the funstools decomposed data table format

from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn5.0d3s1' # the decomposed paramter option title (i.e., directory name)

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option+'_smo','output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_vsm2.0_ssm4.0_v06_cube_{}_g+_fit_fin_finalized.dat'.format(dc_option)
gpy2fitres = 'SerB_C18O_vsm2.0_ssm4.0_v06_cube_{}_g+_fit_fin_finalized_fitres.dat'.format(dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
#        print(j)

    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)

xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, snr, dv, int_tot], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'))
fr.write(os.path.join(path_to_decomp_data, gpy2fitres), format='ascii', overwrite='True')




############# smoothed data 8.25~ GPY+ test
## gpy finalized data  to the funstools decomposed data table format


from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn5.0d3s1' # the decomposed paramter option title (i.e., directory name)
dc_appendix = '_sv2521'
vsm= '5.0'
ssm= '2.0'

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/crop'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option+dc_appendix,'output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_vsm{}_ssm{}_v06_cube_cut_{}_g+_fit_fin_sf-p2_finalized.dat'.format(vsm, ssm, dc_option)
gpy2fitres = 'SerB_C18O_vsm{}_ssm{}_v06_cube_cut_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(vsm, ssm, dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
#        print(j)

    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)

xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, snr, dv, int_tot], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'))
fr.write(os.path.join(path_to_decomp_data, gpy2fitres), format='ascii', overwrite='True')


############# smoothed data 9.04~ GPY+ rectangular map
## gpy finalized data  to the funstools decomposed data table format
# sv2511

from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn3.0d2s1' # the decomposed paramter option title (i.e., directory name)
dc_appendix = '_sv2511'
vsm= '5.0'
ssm= '2.0'

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/rect'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option+dc_appendix,'output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized.dat'.format(vsm, ssm, dc_option)
gpy2fitres = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(vsm, ssm, dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
#        print(j)

    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)

xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, snr, dv, int_tot], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'))
fr.write(os.path.join(path_to_decomp_data, gpy2fitres), format='ascii', overwrite='True')




############# smoothed data 9.04~ GPY+ rectangular map
## gpy finalized data  to the funstools decomposed data table format
# sv3511

from astropy.io.ascii import read, write
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import os

dc_option = 'sn5.0d2s1' # the decomposed paramter option title (i.e., directory name)
dc_appendix = '_sv3511'
vsm= '5.0'
ssm= '3.0'

#gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/ClusterA/'
gpy_dir = '/Users/khkim/Dropbox/Work5/FUNS-Serp/GaussPyPlus/Main/rect'
path_to_decomp_data = os.path.join(
    gpy_dir,dc_option+dc_appendix,'output', 'gpy_decomposed')
gpy_final = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized.dat'.format(vsm, ssm, dc_option)
gpy2fitres = 'SerB_C18O_vsm{}_ssm{}_v06_cube_rect_{}_g+_fit_fin_sf-p2_finalized_fitres.dat'.format(vsm, ssm, dc_option)
#fr = Table(names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'), dtype=('i4', 'i4', 'i4', 'i4', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8'))

dc_data = read(os.path.join(path_to_decomp_data, gpy_final))
cn = []

a = iter(range(len(dc_data)))
for i in a:
#    print('initial i =', i)
#    print(dc_data['ncomps'][i])
    for j in range(dc_data['ncomps'][i]):
        cn.append(j)
#        print(j)

    if dc_data['ncomps'][i] != 1:
        for j in range(dc_data['ncomps'][i]-1):
            next(a)

xpos = np.array(dc_data['x_pos'])
ypos = np.array(dc_data['y_pos'])
ncomps = np.array(dc_data['ncomps'])
amp = np.array(dc_data['amp'])
VLSR = np.array(dc_data['VLSR'])
vel_disp = np.array(dc_data['vel_disp'])
snr = np.array(dc_data['amp']/dc_data['rms'])
dv = np.array(dc_data['vel_disp']*np.sqrt(8.*np.log(2.)))
int_tot = np.array(dc_data['int_tot'])

fr = Table([xpos, ypos, ncomps, cn, amp, VLSR, vel_disp, snr, dv, int_tot], names=('rp', 'dp', 'tn', 'cn', 'tp', 'vp', 'sd', 'snr', 'dv', 'area'))
fr.write(os.path.join(path_to_decomp_data, gpy2fitres), format='ascii', overwrite='True')

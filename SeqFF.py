"'Sequencial FoF (or combined with FIVE)''"
"FOF with maxd <= 2 pixels"
# 1/3/2021

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io.ascii import read, write
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.wcs import WCS
from astropy.table import Table
from mpl_toolkits.mplot3d import axes3d

from mpl_toolkits.mplot3d.axes3d import Axes3D

from astropy.modeling.models import Gaussian1D
from spectral_cube import SpectralCube
from astropy import units as u
import os
from run_funstools_Serp import fof_seq

dc_dir = os.path.join('/Users/khkim/Dropbox/Work5/FUNS-Serp/decompose', 'dc_refit', '4_121921')


########## Start of Step 1 ##############
# read a fof table : snr=2.5, delv=0.12 for the initial selection
dist_v0 = 0.2
dist_rd0 = 2
sigma_cut = 2.5
FoF_dir = dc_dir+'/fof-delv/fof-'+str(dist_v0)+'/'
fof_file = FoF_dir+'FoF_result_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)

new_ff_dir = dc_dir+'/ff_sequence/'

#fn numbers to save separately: independent(?) filaments
# snr=2.5, delv=0.2
selected_fn = [305, 286, 513, 309, 334, 603, 464, 465, 67, 213, 55]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected = ff_table[iselmask]

ff_selected['step'] = 1
ff_selected['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected['fn']==selected_fn[i])
    ff_selected['finfn'][ifmask] = i+1

write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step1.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step1.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step1.dat', overwrite=True)

########## End of Step 1 ##############


########## Start of Step 2 ##############
# read a fof table : snr=2.5, delv=0.1 for the step 2 selection
dist_v0 = 0.16
dist_rd0 = 2
sigma_cut = 2.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step1_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 2 : fn numbers to save separately: independent(?) filaments
# snr=2.5, delv=0.16
selected_fn = [422, 111, 542]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step2.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step2.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step2.dat', overwrite=True)

########## End of Step 2 ##############


########## Start of Step 3 ##############
# read a fof table : snr=2.5, delv=0.08 for the step 3 selection
dist_v0 = 0.1
dist_rd0 = 2
sigma_cut = 2.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step2_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 3 : fn numbers to save separately: independent(?) filaments
# snr=2.5, delv=0.10
selected_fn = [436, 410, 541, 479, 642, 707, 31, 287] #1번 안 떼어냄

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 3
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step3.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step3.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step3.dat', overwrite=True)

########## End of Step 3 ##############


########## Start of Step 4 ##############
# read a fof table : snr=2.5, delv=0.08 for the step 4 selection
dist_v0 = 0.08
dist_rd0 = 2
sigma_cut = 2.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step3_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 4 : fn numbers to save separately: independent(?) filaments
# snr=2.5, delv=0.08
selected_fn = [1, 265, 171, 343, 268, 245, 177, 188, 109]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 4
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step4.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step4.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step4.dat', overwrite=True)



########## Start of Step 5 ##############
# read a fof table : snr=2.5, delv=0.07 for the step 5 selection
dist_v0 = 0.07
dist_rd0 = 2
sigma_cut = 2.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step4_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 4 : fn numbers to save separately: independent(?) filaments
# snr=2.5, delv=0.08
selected_fn = [231, 514, 278]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 5
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step5.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step5.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step5.dat', overwrite=True)




########## Start of Step 6 ##############
# read a fof table : snr=3.0, delv=0.08 for the step 6 selection
dist_v0 = 0.08
dist_rd0 = 2
sigma_cut = 3.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step5_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 4 : fn numbers to save separately: independent(?) filaments
# snr=2.5, delv=0.08
selected_fn = [569, 187, 70, 355, 21, 477, 42, 528]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 6
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step6.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step6.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step6.dat', overwrite=True)




########## Start of Step 7 ##############
# read a fof table : snr=3.0, delv=0.06 for the step 7 selection
dist_v0 = 0.06
dist_rd0 = 2
sigma_cut = 3.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step6_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 7 : fn numbers to save separately: independent(?) filaments
# snr=3.0, delv=0.06
selected_fn = [625, 215, 581, 397, 94, 46, 113, 73, 452]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 7
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step7.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step7.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step7.dat', overwrite=True)



########## Start of Step 8 ##############
# read a fof table : snr=3.0, delv=0.08 again for the step 8 selection
dist_v0 = 0.08
dist_rd0 = 2
sigma_cut = 3.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step7_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 7 : fn numbers to save separately: independent(?) filaments
# snr=3.0, delv=0.06
selected_fn = [66]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 8
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step8.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step8.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step8.dat', overwrite=True)



########## Start of Step 9 ##############
# read a fof table : snr=3.0, delv=0.08 again for the step 8 selection
dist_v0 = 0.05
dist_rd0 = 2
sigma_cut = 3.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step8_remain_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 9 : fn numbers to save separately: independent(?) filaments
# snr=3.0, delv=0.05
selected_fn = [87, 120, 19, 71, 171, 321, 50]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = 9
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step9.dat', overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step9.dat', overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step9.dat', overwrite=True)



########## Start of Step 10 ##############
# read a fof table : snr=3.5, delv=0.12 again for the step 9selection
step1=9
step2=10
dist_v0 = 0.12
dist_rd0 = 2
sigma_cut = 3.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 10 : fn numbers to save separately: independent(?) filaments
# snr=3.5, delv=0.12
selected_fn = [215, 55]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)



########## Start of Step 11 ##############
# read a fof table : snr=3.5, delv=0.08 again for the step 9selection
step1=10
step2=11
dist_v0 = 0.08
dist_rd0 = 2
sigma_cut = 3.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 11 : fn numbers to save separately: independent(?) filaments
# snr=3.5, delv=0.08
selected_fn = [88]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)



########## Start of Step 12 ##############
# read a fof table : snr=3.5, delv=0.06 again for the step 9selection
step1=11
step2=12
dist_v0 = 0.06
dist_rd0 = 2
sigma_cut = 3.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 12 : fn numbers to save separately: independent(?) filaments
# snr=3.5, delv=0.06
selected_fn = [24]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)





########## Start of Step 13 ##############
# read a fof table : snr=3.5, delv=0.06 again for the step 9selection
step1=12
step2=13
dist_v0 = 0.12
dist_rd0 = 2
sigma_cut = 4.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 12 : fn numbers to save separately: independent(?) filaments
# snr=3.5, delv=0.06
selected_fn = [120, 19, 130, 205, 89]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)



########## Start of Step 14 ##############
# read a fof table : snr=4.0, delv=0.1 again for the step 14 selection
step1=13
step2=14
dist_v0 = 0.1
dist_rd0 = 2
sigma_cut = 4.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 14 : fn numbers to save separately: independent(?) filaments
# snr=4.0, delv=0.1
selected_fn = [118]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)


########## Start of Step 15 ##############
# read a fof table : snr=4.0, delv=0.1 again for the step 14 selection
step1=14
step2=15
dist_v0 = 0.07
dist_rd0 = 2
sigma_cut = 4.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 14 : fn numbers to save separately: independent(?) filaments
# snr=4.0, delv=0.1
selected_fn = [221]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)


########## Start of Step 16 ##############
# read a fof table : snr=4.0, delv=0.1 again for the step 14 selection
step1=15
step2=16
dist_v0 = 0.08
dist_rd0 = 2
sigma_cut = 4.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 14 : fn numbers to save separately: independent(?) filaments
# snr=4.0, delv=0.1
selected_fn = [216, 166]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)


########## Start of Step 17 ##############
# read a fof table : snr=4.5, delv=0.07 again for the step 17 selection
step1=16
step2=17
dist_v0 = 0.07
dist_rd0 = 2
sigma_cut = 4.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 17 : fn numbers to save separately: independent(?) filaments
# snr=4.5, delv=0.07
selected_fn = [209, 35]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)



########## Start of Step 18 ##############
# read a fof table : snr=4.5, delv=0.07 again for the step 17 selection
step1=17
step2=18
dist_v0 = 0.06
dist_rd0 = 2
sigma_cut = 4.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 17 : fn numbers to save separately: independent(?) filaments
# snr=4.5, delv=0.07
selected_fn = [24, 42, 44]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)






########## Start of Step 19 ##############
# read a fof table : snr=4.5, delv=0.05 again for the step 17 selection
step1=18
step2=19
dist_v0 = 0.05
dist_rd0 = 2
sigma_cut = 4.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 17 : fn numbers to save separately: independent(?) filaments
# snr=4.5, delv=0.07
selected_fn = [31, 39]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)



########## Start of Step 20 ##############
# read a fof table : snr=4.5, delv=0.05 again for the step 17 selection
step1=19
step2=20
dist_v0 = 0.06
dist_rd0 = 2
sigma_cut = 5.0
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 17 : fn numbers to save separately: independent(?) filaments
# snr=4.5, delv=0.07
selected_fn = [67]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)



########## Start of Step 21 ##############
# read a fof table : snr=4.5, delv=0.05 again for the step 17 selection
step1=20
step2=21
dist_v0 = 0.2
dist_rd0 = 2
sigma_cut = 5.5
FoF_dir = dc_dir+'/ff_sequence/'
fof_file = FoF_dir+'process/FoF_result_C18O_SeqFF_step{}_remain_'.format(step1)+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'.dat'
ff_table = read(fof_file)
SeqFF_file = FoF_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat'
ff_selected = read(SeqFF_file)
maxfn = np.max(ff_selected['finfn'])

new_ff_dir = dc_dir+'/ff_sequence/'

#Step 17 : fn numbers to save separately: independent(?) filaments
# snr=4.5, delv=0.07
selected_fn = [1, 2, 15]

#subtable with the selected fn
selmask=[]
for i in range(np.size(selected_fn)):
    selmask.append(np.where(ff_table['fn']==selected_fn[i])[0])
iselmask = np.concatenate(selmask).ravel()
ff_selected_now = ff_table[iselmask]

ff_selected_now['step'] = step2
ff_selected_now['finfn'] = 0

for i in range(np.size(selected_fn)):
    ifmask = np.where(ff_selected_now['fn']==selected_fn[i])
    ff_selected_now['finfn'][ifmask] = maxfn+i+1

for i in range(np.size(ff_selected_now)):
    ff_selected.add_row(ff_selected_now[i])

write(ff_selected_now, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4_step{}.dat'.format(step2), overwrite=True)
write(ff_selected, new_ff_dir+'ff_selected/'+'SeqFF_selected_C18O_dc_comb_rev4.dat', overwrite=True)


#the remains for the later fof
ff_keep = ff_table.copy()
ff_keep.remove_rows(iselmask)
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_step{}.dat'.format(step2), overwrite=True)

ff_keep['fn'] = 0
write(ff_keep, new_ff_dir+'ff_keep/'+'SeqFF_remain_C18O_dc_comb_rev4_'+str(float(sigma_cut))+'sig_v'+str(dist_v0)+'_d'+str(dist_rd0)+'_fn0_step{}.dat'.format(step2), overwrite=True)

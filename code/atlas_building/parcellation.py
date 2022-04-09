#%matplotlib qt5
#import matplotlib.pyplot as plt
from nilearn import plotting
from nilearn.image import mean_img, index_img, concat_imgs, resample_img
from nilearn.input_data import NiftiMasker
from nilearn.decomposition import CanICA
import numpy as np
import os

from nilearn.regions import Parcellations
import nibabel as nib

source ='source/'
target = 'output/'
exclude = ['empty','SEZL','SEZR','WB','other']
include = ['nii']

# Parameters
mp_comp    = 80
mp_smooth  =  5
mp_detrend = False
mp_std     = False
mp_jobs    = 4
mp_method  = 'ica'
mp_label   = 'average_dff_SYH_2states_fullMask'
downsample = 2
mp_affine  = np.eye(4)*(downsample,downsample,1,1)


out_name = ''.join([mp_method, '_', mp_label,
                                    '_comps', str(mp_comp),
                                    '_smooth', str(mp_smooth),
                                    '_detrend', str(mp_detrend),
                                    '_std', str(mp_std),
                                    '_resampled', str(downsample)])

print('Will be saving to ' + target + out_name + '.nii')

files_list = list()
for path, subdirs, files in os.walk(source):
    for name in files:
        fullpath = os.path.join(path, name)
        files_list.append(fullpath)

files_list = [f for f in files_list if '.nii' in f and any(inc in f for inc in include)]
files_list = [f for f in files_list if '.nii' in f and not any(ex in f for ex in exclude)]
print(files_list)

SEZmask = 'SEZ_template_2.0_mask_int8.nii'

import time

canica = CanICA(n_components=mp_comp, smoothing_fwhm=mp_smooth,
                memory="/tmp/nilearn_cache", memory_level=2, n_jobs=mp_jobs,
                threshold=None, verbose=100, random_state=None, standardize=mp_std, detrend=mp_detrend, mask=SEZmask, target_affine=mp_affine)

start = time.time()
canica.fit(files_list)
print("CanICA " + str(mp_comp) + "clusters: %.2fs" % (time.time() - start))
components_img = canica.components_img_
nib.save(nib.Nifti1Image(np.float32(components_img.get_fdata()), np.eye(4)), target + out_name + '.nii')

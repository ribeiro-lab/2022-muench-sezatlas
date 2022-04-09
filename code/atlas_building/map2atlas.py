# run region extractor, assign each voxel to 1 region, cleanup and erode

import matplotlib.pyplot as plt
from nilearn import plotting
from nilearn.image import mean_img, index_img, concat_imgs, resample_img
from nilearn.regions import RegionExtractor
import numpy as np
import os
import nibabel as nib

header = nib.load('SEZ_template_2.0.nii.gz').header
source_dir = 'source/'
output     = 'output/'

maps = [m for m in os.listdir(source_dir) if m.endswith('nii')]


def disk_3d(radius):
    a = radius
    test = sp.zeros((a*2+1,a*2+1,a*2+1), dtype=sp.int16)
    test[:,:,a] = disk(a)
    return test

for map in maps: #map = maps[0]
    print('\n\nPROCESSING: ' + map)
    input = source_dir + map

    ######################
    ## REGION EXTRACTOR ##
    ######################

    # Parameters
    mp_thres     = 1
    mp_smooth    = 8
    mp_size      = 1350
    mp_extractor = 'connected_components' #'local_regions'#
    mp_strategy  = 'ratio_n_voxels'
    mp_label     = map[:-4] + '_gauss2200'

    print('     Will be saving to: ' + mp_label + '.nii')

    extractor = RegionExtractor(input, threshold=mp_thres,
                                smoothing_fwhm=mp_smooth,
                                thresholding_strategy=mp_strategy,
                                extractor=mp_extractor,
                                standardize=False,
                                min_region_size=mp_size,
                                verbose=10)

    extractor.fit()
    regions_extracted_img = extractor.regions_img_
    regions_index = extractor.index_
    n_regions_extracted = regions_extracted_img.shape[-1]

    data = np.float32(regions_extracted_img.get_fdata())
    # gauss filter
    print('applying gaussian filter on components')
    from scipy.ndimage import gaussian_filter
    data = gaussian_filter(data, sigma=(2,2,0,0))

    nifti_img = nib.Nifti1Image(data, np.eye(4), header = header)
    nifti_img.header['xyzt_units']
    nifti_img.header['pixdim'] = header['pixdim']

    nib.save(nifti_img, output + 'regions_' + mp_label + '_regions_' + str(n_regions_extracted) + '.nii')
    print("     saved regions...")


    ##################
    ## ASSIGN LABEL ##
    ##################

    import copy

    data = nifti_img.get_fdata()

    def assign_label(regions):
        '''Assignes a voxel to the region/label with the strongest contribution'''
        max_indx = np.where(regions == max(regions))[0]
        if len(max_indx) > 1:
            raise ValueError('More than one maximal value!')
        regions.fill(0)
        regions[max_indx[0]] = max_indx[0] + 1
        return regions

    # reshape data
    x,y,z,c = data.shape
    data_2D = data.reshape(x*y*z,c)
    data_2D.shape

    # find voxels that contribute to several components
    sum = copy.copy(data_2D)
    indx = sum > 0
    sum[indx] = 1
    sum.shape
    sum = np.sum(sum, axis=1)
    sum.shape
    indx_l1 = np.where(sum > 1)[0]
    indx_1 = np.where(sum == 1)[0]
    len_l1 = len(sum[indx_l1])
    len_1 = len(sum[indx_1])
    perc = round(float(len_l1) / float(len_1) * 100,1)
    print('     ', str(perc) + '% of ' + str(len_1) + ' active voxels contribute to more than 1 component.')

    # assign label id
    for i in indx_l1:
        data_2D[i,:] = assign_label(data_2D[i,:])

    for i in indx_1:
        data_2D[i,:] = assign_label(data_2D[i,:])

    # flatten array
    data_flat = np.amax(data_2D, axis=1)
    new_data = data_flat.reshape(x,y,z)
    labels = nib.Nifti1Image(np.int16(new_data), np.eye(4), header=nifti_img.header)
    nib.save(labels, output + mp_label +'_atlas_raw.nii')
    print('     saved raw atlas...')


    #############
    ## CLEANUP ##
    #############

    import scipy as sp
    import scipy.ndimage as ndimage
    from skimage.morphology import disk, ball
    header = labels.header
    data = labels.get_fdata()

    #remove regions smaller 5x5x5 (125) pixels
    unique, counts = np.unique(data, return_counts = True)
    np.where(counts <= 125)
    low = unique[np.where(counts <= 125)]
    for i in low:
        data[np.where(data == i)] = 0

    print("regions < 125: " + str(low))

    values = list(np.unique(data))
    values.remove(0)

    opened = sp.zeros(data.shape, dtype=sp.int16)
    for value in values:
        subset = sp.zeros(data.shape, dtype=sp.int16)
        subset[sp.where(data == value)] = value
        subset = sp.ndimage.binary_opening(subset, structure = ball(1.5))
        opened[sp.where(subset == True)] = value

    filled = sp.zeros(data.shape, dtype=sp.int16)
    for value in values:
        subset = sp.zeros(opened.shape, dtype=sp.int16)
        subset[sp.where(opened == value)] = value
        subset = ndimage.binary_fill_holes(subset, structure = ball(3)) #sp.ones((5,5,5)))
        filled[sp.where(subset == True)] = value


    closed = sp.zeros(data.shape, dtype=sp.int16)
    for value in values:
        subset = sp.zeros(filled.shape, dtype=sp.int16)
        subset[sp.where(filled == value)] = 1
        #plt.imshow(subset[:,:,24].T)
        #subset = sp.ndimage.binary_closing(subset, structure = sp.ones((3,3,3)))
        subset = ndimage.binary_closing(subset, structure = ball(1.5))
        closed[sp.where(subset == True)] = value

    closed2 = sp.zeros(data.shape, dtype=sp.int16)
    for value in values:
        subset = sp.zeros(closed.shape, dtype=sp.int16)
        subset[sp.where(closed == value)] = 1
        #plt.imshow(subset[:,:,24].T)
        #subset = sp.ndimage.binary_closing(subset, structure = sp.ones((3,3,3)))
        subset = ndimage.binary_closing(subset, structure = disk_3d(4))
        closed2[sp.where(subset == True)] = value

    filled2 = sp.zeros(data.shape, dtype=sp.int16)
    for value in values:
        subset = sp.zeros(closed2.shape, dtype=sp.int16)
        subset[sp.where(closed2 == value)] = value
        subset = ndimage.binary_fill_holes(subset, structure = disk_3d(4)) #sp.ones((5,5,5)))
        filled2[sp.where(subset == True)] = value

    eroded = sp.zeros(data.shape, dtype=sp.int16)
    for value in values:
        subset = sp.zeros(filled2.shape, dtype=sp.int16)
        subset[sp.where(filled2 == value)] = value
        subset = sp.ndimage.binary_erosion(subset, structure = disk_3d(1))
        eroded[sp.where(subset == True)] = value

    final = eroded

    unique, counts = sp.unique(final, return_counts = True)
    sp.where(counts <= 125)
    low = unique[sp.where(counts <= 125)]
    for i in low:
        final[sp.where(final == i)] = 0

    print("regions < 125: " + str(low))

    unique, counts = sp.unique(final, return_counts = True)
    j = 0
    for i in unique:
        final[sp.where(final == i)] = j
        j = j+1


    atlas = nib.Nifti1Image(sp.int16(final), sp.eye(4), header=header)
    nib.save(atlas, output + mp_label + '_atlas' +'.nii')
    print('     saved final atlas...\n\n\n')

import numpy as np
import os

from astropy import units as u
from astropy.io import fits
from astropy.wcs import WCS
import gausspyplus.plotting as gplt 
from tqdm import tqdm
import math
import itertools

import pickle
import random
import sys

import matplotlib
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from astropy import units as u
from gausspyplus.utils.spectral_cube_functions import correct_header, change_header, update_header

def pickle_load(path_to_file, binary=True, encoding='latin1'):
    import pickle
    read = 'r'
    if binary:
        read = 'rb'
    with open(path_to_file, read) as pickled_file:
        pickled_data = pickle.load(pickled_file, encoding=encoding)
    return pickled_data

def getgauss( pixel_range = None):
    
    vel_unit=u.km/u.s
    cols=5
    rowsize=7.75
    rowbreak=50
    dpi=50
    n_spectra=None
    suffix=''
    subcube=False
    list_indices=None
    gaussians=True
    residual=True
    
    random_seed=111
    
    channels, fig_channels = (data['x_values'] for _ in range(2))
    n_channels = len(channels)
    
    header = gplt.correct_header(data['header'])
    
    
    list_indices, n_spectra, grid_layout = gplt.get_list_indices(data, subcube=subcube, pixel_range=pixel_range, list_indices=list_indices, n_spectra=n_spectra, 
                                                                 random_seed=random_seed)
    cols, rows, rowbreak, colsize, multiple_pdfs = gplt.get_figure_params(n_channels, n_spectra, cols, rowsize, rowbreak, grid_layout, subcube=subcube)
    
    pbar = tqdm(total=n_spectra)
    
    for i, idx in enumerate(list_indices):
        pbar.update(1)
        yi, xi = data['location'][idx]
        index_data = data['index'][idx]
        
        if rowbreak is not None:
            k = int(i / (rowbreak*cols))
            if (k + 1)*rowbreak > rows:
                rows_in_figure = rows - k*rowbreak

            else:
                rows_in_figure = rowbreak

        y = data['data_list'][idx]
        rms = data['error'][idx][0]
        
        fit_fwhms = decomp['fwhms_fit'][idx]
        fit_means = decomp['means_fit'][idx]
        fit_amps = decomp['amplitudes_fit'][idx]
        
        gauss = gplt.combined_gaussian(fit_amps, fit_fwhms, fit_means, channels)
        
        
    return(gauss)

def getdisk(xm,ym, pixel_range = None):
    
    n_spectra=None
    subcube=False
    list_indices=None
    gaussians=True
    
    random_seed=111
    
    channels, fig_channels = (data['x_values'] for _ in range(2))
    n_channels = len(channels)
    
    
    
    list_indices, n_spectra, grid_layout = gplt.get_list_indices(data, subcube=subcube, pixel_range=pixel_range, list_indices=list_indices, 
                                                                 n_spectra=n_spectra, random_seed=random_seed)

    
    for i, idx in enumerate(list_indices):

        yi, xi = data['location'][idx]
        index_data = data['index'][idx]

        y = data['data_list'][idx]
        rms = data['error'][idx][0]
        
        fit_fwhms = decomp['fwhms_fit'][idx]
        fit_means = decomp['means_fit'][idx]
        fit_amps = decomp['amplitudes_fit'][idx]

        gauss = np.zeros(189)
        min_diff = 189
        minj = 0
        ncomps = len(fit_amps)
        
        max_chan = 0
        max_flux = 0

        for i in range(len(bb_image[:,yi,xi])):
            if bb_image[i,yi,xi] > max_flux :
                max_chan = i
                max_flux = bb_image[i,yi,xi]
        
        for j in range(ncomps):
            if ((np.abs(max_chan-fit_means[j])) < min_diff) and max_flux > 0:
                min_diff = np.abs(max_chan-fit_means[j])
                minj = j
                gauss = gplt.gaussian(fit_amps[minj], fit_fwhms[minj], fit_means[minj], channels)
            
        
    return(gauss)

def getcomps(data,decomp,xm,ym, model, cube, mom0, vsys, pixel_range = None):
    
    n_spectra=None
    subcube=False
    list_indices=None
    gaussians=True
    random_seed=111
    dv = 8544.9290-8539.4159

    
    channels, fig_channels = (data['x_values'] for _ in range(2))
    n_channels = len(channels)
    
    
    list_indices, n_spectra, grid_layout = gplt.get_list_indices(data, subcube=subcube, pixel_range=pixel_range, list_indices=list_indices, 
                                                                 n_spectra=n_spectra, random_seed=random_seed)


    for i, idx in enumerate(list_indices):

        yi, xi = data['location'][idx]
        index_data = data['index'][idx]

        y = data['data_list'][idx]
        rms = data['error'][idx][0]
        
        fit_fwhms = decomp['fwhms_fit'][idx]
        fit_means = decomp['means_fit'][idx]
        fit_amps = decomp['amplitudes_fit'][idx]

        gaussdisk = np.zeros(189)
        gausslead = np.zeros(189)
        gaussfollow = np.zeros(189)
        gaussanom = np.zeros(189)
        anommom = 0
        min_diff = 189
        minj = 0
        ncomps = len(fit_amps)
        
        max_chan = 0
        max_flux = 0
        vel = 0

        means = np.empty(ncomps)
        tag = np.empty(ncomps)
        diskcomp = 0
        
        for i in range(len(model[:,yi,xi])):
            if model[i,yi,xi] > max_flux  :
                max_chan = i
                max_flux = model[i,yi,xi]
        
        for j in range(ncomps):
            if ((np.abs(max_chan-fit_means[j])) < min_diff) and max_flux > 0:
                min_diff = np.abs(max_chan-fit_means[j])
                minj = j
                gaussdisk = gplt.gaussian(fit_amps[minj], fit_fwhms[minj], fit_means[minj], channels)
                diskcomp = fit_means[minj]
                
                
        direction = (8544.9290 - diskcomp*dv) - vsys

        if direction > 0 :
            for j in range(ncomps):
                comp = fit_means[j]-diskcomp
                if comp > 0 and max_flux > 0:
                    gausslead = gplt.gaussian(fit_amps[j], fit_fwhms[j], fit_means[j], channels)
                elif comp < 0 and max_flux > 0:
                    gaussfollow = gplt.gaussian(fit_amps[j], fit_fwhms[j], fit_means[j], channels)
                
        else:
            for j in range(ncomps):
                comp = fit_means[j]-diskcomp
                if comp < 0 and max_flux > 0:
                    gausslead = gplt.gaussian(fit_amps[j], fit_fwhms[j], fit_means[j], channels)
                elif comp > 0 and max_flux > 0:
                    gaussfollow = gplt.gaussian(fit_amps[j], fit_fwhms[j], fit_means[j], channels)
                
        if ncomps == 0:
            gaussanom = cube[:,yi,xi]
            anommom = mom0[yi,xi]
            
        
    return(gaussdisk, gausslead, gaussfollow, gaussanom, anommom)
    
    
def make_images(data,
                decomp,
                cube_image,
                cube_mom0,
                v,
                cube,
                image):
    
    cube_copy = np.zeros_like(cube_image)
    cube_lead = np.zeros_like(cube_image)
    cube_follow = np.zeros_like(cube_image)
    cube_anom = np.zeros_like(cube_image)
    anom_mom0 = np.zeros_like(cube_mom0)


    s = len(cube_copy[:,0,0])
    x = len(cube_copy[0,:,0])
    y = len(cube_copy[0,0,:])
    k = 0

    for i in range(0,x):
        for j in range(0,y):
        
            pixel_range =   {'x': [j,j], 'y': [i,i]}
            gaussdisk, gausslead, gaussfollow, gaussanom,anommom = getcomps(data,decomp,j, i, image, cube, cube_mom0, v, pixel_range = pixel_range)
            cube_copy[:,i,j] = gaussdisk
            cube_lead[:,i,j] = gausslead
            cube_follow[:,i,j] = gaussfollow
            cube_anom[:,i,j] = gaussanom
            anom_mom0[i,j] = anommom
            
    return cube_copy, cube_lead, cube_follow, cube_anom, anom_mom0
        
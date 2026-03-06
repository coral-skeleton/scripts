import numpy as np
import matplotlib.pyplot as plt
import math
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

expected_lines = {
    'OII_3726': 3726,
    'OII_3729': 3729,
    'NeIII_3869': 3869,
    'CaH': 3934,
    'CaK': 3968,
    'NeIII_3967': 3967,
    'Hδ': 4101,
    'Hγ': 4340,
    'Hb': 4861.333,
    'OIII_4959': 4958.911,
    'OIII_5007': 5006.843,
    'Mgb':5175,
    'NI_5198': 5198,
    'NI_5200': 5200,
    'Fe_5270': 5270,
    'Fe_5335': 5335,
    'HeI_5876': 5876,
    'NaD_5890': 5890,
    'NaD_5896': 5896,
    'OI_6300': 6300,
    'OI_6363': 6363,
    'SIII': 6312, 
    'NII_6548': 6548.05,
    'Ha': 6562.819,
    'NII_6583': 6583.46,
    'HeI_6678': 6678,
    'SII_6716': 6716.44,
    'SII_6731': 6730.81,
    'HeI_7065': 7065,
}

def chipfill(data, refspec): #interpolating over the chip gap in 2d rss image, data and refspec are numpy arrays containing images
    
    chip_fill= np.copy(data)
    x,y = data.shape
    l = np.arrange(0,y,num=y)
    for i in range(0,x):        
        m, c = np.polyfit(real, data[i,:], 1)
        y = m*l + c
        for j in range(0,len(refspec[0,10:y+10])):
            if refspec[i,(j+10)] <= 1:
                chip_fill[i,(j+10)] = y[(j+10)]
                chip_fill[i,(j+10)+1] = y[(j+10)+1]
                chip_fill[i,(j+10)-1] = y[(j+10)-1]
                chip_fill[i,(j+10)+2] = y[(j+10)+2]
                chip_fill[i,(j+10)-2] = y[(j+10)-2]
                
    return chip_fill


def flatred(data,flat): #flat fielding for rss data, data and flat are numpy arrays containing images
    
    spatial_profile = np.median(flat, axis=1)
    smooth_profile = gaussian_filter(spatial_profile, sigma=20)
    illum_profile = smooth_profile / np.median(smooth_profile)
    illum_flat = illum_profile[:,None]

    illum_flat /= np.mean(illum_flat)
    flat_red = data / illum_flat
    
    return flat_red

def stackimages(images = []): #for stacking a set of images, input is a list of numpy arrays containing images, eg [data1,data2,data3]

    stack = np.stack(images, axis=0)  

    stacked = np.sum(stack, axis=0)

    sigma = np.std(stack, axis=0, ddof=1) 
    
    return stacked, sigma
    
def simpleskysubifu(data, skyfib = [4,22,30,31,40,52,64,77,91,104,117,132,155,167,182,194,237,239,242,250,275,291,300,309]):
#simple skysub for ifu
    sky_fibers = np.array(skyfib) 
    sky_spectra = data[sky_fibers - 1]
    clipped = sigma_clip(sky_spectra, axis=0, sigma=3)
    master_sky = np.nanmedian(clipped, axis=0)
    skysub = data - master_sky[None, :]
    
    return skysub

def stdskysubifu(data, skyfib = [4,22,30,31,40,52,64,77,91,104,117,132,155,167,182,194,237,239,242,250,275,291,300,309]):
#standard skysub routine that mirrors iraf sky subtraction for ifu
    n_fibers, n_wave = data.shape
    fiber_numbers = np.arange(n_fibers)
    sky_loc = np.array(skyfib)
    continuum_sub, continuum_model = contsub(data)
    sky_model = np.zeros_like(data)

    for w in range(n_wave):

        y_sky = continuum_sub[sky_loc-1, w]
        x_sky = sky_loc-1
        coeffs = np.polyfit(x_sky, y_sky, deg=4)
        poly = np.poly1d(coeffs)

        sky_model[:, w] = poly(fiber_numbers)

    line_subtracted = continuum_sub - sky_model

    skysub = line_subtracted + continuum_model
    
    return skysub, sky_model


def getskyrows(data): #autodetection of longslit sky rows, doesnt care about object trace
    
    nedge=30 
    sigma=3
    
    n_spatial, _ = data.shape

    edge_rows = np.concatenate([
        np.arange(0, nedge),
        np.arange(n_spatial - nedge, n_spatial)
    ])

    row_flux = np.nanmedian(data, axis=1)

    sky_flux = row_flux[edge_rows]

    med = np.median(sky_flux)
    std = np.std(sky_flux)

    good = np.abs(sky_flux - med) < sigma*std

    skyrows = edge_rows[good]
    
    return skyrows


def getskyrows_auto(data, sigma=3, smooth=3, padding=5): #automatic detection of longslit sky rows, keeping object trace in mind

    n_spatial, n_wave = data.shape

    spatial_profile = np.nanmedian(data, axis=1)

    smooth_profile = gaussian_filter1d(spatial_profile, smooth)

    med = np.median(smooth_profile)
    std = np.std(smooth_profile)

    object_rows = np.where(smooth_profile > med + sigma*std)[0]

    if len(object_rows) == 0:
        return np.arange(n_spatial)

    obj_min = max(object_rows.min() - padding, 0)
    obj_max = min(object_rows.max() + padding, n_spatial)

    skyrows = np.concatenate([
        np.arange(0, obj_min),
        np.arange(obj_max, n_spatial)
    ])

    return skyrows

def simpleskysubls(data): #simple skysub for longslit
    
    skyrows = getskyrows(data)
    
    sky_spectra = data[skyrows]

    clipped = sigma_clip(sky_spectra, axis=0, sigma=3)

    master_sky = np.nanmedian(clipped, axis=0)

    skysub = data - master_sky[None, :]
    
    return skysub

def stdskysubls(data): #standard skysub routine that mirrors iraf sky subtraction for longslit
    
    skyrows = getskyrows(data)
    n_spatial, n_wave = data.shape

    spatial_pixels = np.arange(n_spatial)
    sky_loc = np.array(skyrows)

    continuum_sub, continuum_model = contsub(data)

    sky_model = np.zeros_like(data)

    for w in range(n_wave):

        y_sky = continuum_sub[sky_loc, w]
        x_sky = sky_loc

        coeffs = np.polyfit(x_sky, y_sky, deg=4)
        poly = np.poly1d(coeffs)

        sky_model[:, w] = poly(spatial_pixels)

    line_subtracted = continuum_sub - sky_model

    skysub = line_subtracted + continuum_model
    
    return skysub, sky_model

def maskspectrum(data, l, minl, maxl, w=10): #spectral masking, l is redshift corrected spectral range, lmin and lmax is min and max of wavelenght range, only works for wavelenght range between 3500-7100A
    
    cent, continuum = contsub(data)
    
    width = w
    mask = np.ones_like(cent, dtype=bool)

    for label, line_center in expected_lines.items():
        if line_center >= minl and line_center <= maxl:
            mask[:, abs(l - line_center) < width] = False

    masked = cent.copy()
    masked[~mask] = 0
    
    return masked

def contsub(data): #continuum subtraction
    
    x,y = data.shape
    wl = np.rouind(2*(y/3))
    if wl%2 == 0:
        wl += 1
    continuum = savgol_filter(data, window_length=wl, polyorder=4)
    cent = data-continuum
    
    return cent, continuum

def pca(data,l,lmin,lmax): #singular value decomp of data
    
    
    masked = maskspectrum(data,l,lmin,lmax)
    
    U, S, Vt = np.linalg.svd(masked, full_matrices=False)
    
    return U, S, Vt

    
    
def idpcacomps(data,l,lmin,lmax):#plot first components for user identification of which to use
    
    x,y = data.shape
    n = np.round(x/2)
    U,S,Vt = pca(data,l,lmin,lmax)
    cent,continuum = contsub(data)
    fig = plt.figure(figsize = (30,8))
    plt.plot(real,cent[n,:]-0.2,label='data')    


    for k in range(12):
        plt.plot(real, Vt[k]+0.2*k,label='component '+str(k))
    
    plt.legend(fontsize=12, ncol = 13)
    plt.title('Component spectra')
    
    variance_fraction = S**2 / np.sum(S**2)
    
    fig = plt.figure(figsize = (8,8))
    plt.plot(variance_fraction, marker='o')
    plt.yscale('log')
    plt.xlabel("Component")
    plt.ylabel("Fraction of Variance")
    plt.title("Variance of components")

def pcaskysub(data, nsky = [1,2,3,4,5],l,lmin,lmax): #actual pca skysub routine
    
    U,S,Vt = pca(data,l,lmin,lmax)
    sky_model = np.zeros_like(data)

    for k in nsky:
        sky_model += np.outer(U[:,k] * S[k], Vt[k])
    
    sci_sub = data - sky_model
    
    return sky_model, skysub


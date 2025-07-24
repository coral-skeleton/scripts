from astroscrappy import detect_cosmics
import numpy as np
import matplotlib.pyplot as plt
import math
#import pyraf
import astropy.units as u
from astropy.utils.data import download_file
from astropy.io import fits  # We use fits to open the actual data file

from astropy.utils import data
data.conf.remote_timeout = 60

from astropy.wcs import wcs
from astropy.visualization import quantity_support
from scipy.ndimage import uniform_filter1d
from scipy.signal import savgol_filter
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.patches as patches
from matplotlib.transforms import Bbox


def normz(gal, sigma, l, expected_lines, Ha = 6562.819, knownz = False, z = 0):
    
    m, c = np.polyfit(l, gal, 1)
    y2 = m*l + c
    gal2 = np.copy(gal)

    for i in range(len(l)):
        if (l[i] < 5840 and l[i] > 5775) or (l[i] > 6817 and l[i] < 6880):
            gal2[i] = y2[i]
        if (gal[i] <= 1e-18):
            gal2[i] = y2[i]
            if (l[i] > 7500):
                for j in range(2):
                    gal2[i-(1+j)] = y2[i-(1+j)]
            else:
                for j in range(2):
                    gal2[i+(1+j)] = y2[i+(1+j)]
                    gal2[i-(1+j)] = y2[i-(1+j)]
                    

    continuum = savgol_filter(gal2, window_length=3171, polyorder=3)
    galnorm = gal2 - continuum
                
    rms = np.sqrt(np.mean(galnorm**2))

    m, c = np.polyfit(l, sigma, 1)
    sigy = m*l + c

    for i in range(len(l)):
        if  sigma[i] <= 0.0:
            sigma[i] = sigy[i]
            if (l[i] < 6000 and l[i] > 5500) or (l[i] > 6500 and l[i] < 7000):
                for j in range(6):
                    sigma[i+1+j] = sigy[i]
                    sigma[i-(1+j)] = sigy[i]
            if (l[i] > 7500):
                for j in range(6):
                    sigma[i-(1+j)] = sigy[i]

    contsig = savgol_filter(sigma, window_length=3171, polyorder=3)
    signorm = sigma-contsig                
     
        
    if knownz:
        real = l/(z+1)
    else:
        ha = 0
        hi = 0

        for i in range(3172):
            if galnorm[i] > hi:
                hi = galnorm[i]
                ha = l[i]
        z = (ha - Ha)/Ha
        real = l/(z+1)

    rms = np.sqrt(np.mean(galnorm ** 2))
    known_lines = np.array(list(expected_lines.values()))
    exclusion_width = 25  # Å
    window_size = 50 
    squared = galnorm ** 2
    rolling_rms = np.sqrt(uniform_filter1d(squared, size=window_size))

    effective_noise = np.maximum(3*signorm, rolling_rms)

    for i in range(len(real)):
        for j in range(len(known_lines)):
            if (real[i] <= known_lines[j]*1.0005 and real[i] >= known_lines[j] - known_lines[j]*0.0005):
                for k in range(exclusion_width):
                    effective_noise[i+k] = 3*sigma[i+k]
                    effective_noise[i-k] = 3*sigma[i-k]

    
    #effective_noise = 3*signorm
    snr = galnorm/effective_noise
    
    return real, galnorm, snr, z, effective_noise
    
    
def plotpeaks(ax, real, galnorm, line_peaks, linenames, linelist):
    ax.plot(real[line_peaks], galnorm[line_peaks], 'rx')

    # Label each peak
    used_label_boxes = []
    for i in line_peaks:
        x = real[i]
        y = galnorm[i]
    
        # Match to line name or fallback
        label = None
        for j in range(len(linelist)):
            wl = float(linelist[j])
            if abs(x - wl) < wl * 0.0005 and y > 0:
                label = linenames[j]
                break
        if label is None:
            label = f"{x:.2f}"
    
        place_label(ax, x, y, label, used_label_boxes)
        
def is_bbox_overlapping(new_bbox, used_label_boxes):
    """Check if new bounding box overlaps any used ones."""
    for box in used_label_boxes:
        if box.overlaps(new_bbox):
            return True
    return False

# Function to try different label placements
def place_label(ax, x, y, text, used_label_boxes):
    """Place a label avoiding overlaps using arrow annotation."""
    offsets = [
        (0, 2e-18),     # above
        (8, 1e-18),     # up and right
        (-8, 1e-18),    # up and left
        (12, 0),        # far right
        (-12, 0),       # far left
        (0, 4e-18),     # higher
    ]
    
    for dx, dy in offsets:
        ann = ax.annotate(
            text,
            xy=(x, y),
            xytext=(x + dx, y + dy),
            textcoords='data',
            arrowprops=dict(arrowstyle='-', color='gray', lw=0.5),
            fontsize=16,
            ha='center',
            va='bottom',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1.0)
        )
        # Draw it to get its position
        fig = ax.get_figure()
        fig.canvas.draw()
        
        # Get label bounding box in data coords
        renderer = fig.canvas.get_renderer()
        bbox = ann.get_window_extent(renderer=renderer)
        bbox_data = bbox.transformed(ax.transData.inverted())

        if not is_bbox_overlapping(bbox_data, used_label_boxes):
            used_label_boxes.append(bbox_data)
            return

        # If overlap, remove this attempt and try next
        ann.remove()
        
def dogauss(real, galnorm, line_peaks):
    fit_results = []

    for i in line_peaks:
        x0 = real[i]
        y0 = galnorm[i]
    
        # Define fitting window: ±10 Å around the peak
        window = 10
        mask = (real > x0 - window) & (real < x0 + window)
        x_fit = real[mask]
        y_fit = galnorm[mask]

        # Initial guess: amp, center, sigma, offset
        amp_guess = y0
        cen_guess = x0
        sigma_guess = 3.0  # Å — typical for SALT
        offset_guess = np.median(y_fit[:3])

        try:
            popt, pcov = curve_fit(
                gaussian,
                x_fit,
                y_fit,
                p0=[amp_guess, cen_guess, sigma_guess, offset_guess]
            )

            amp, cen, sigma, offset = popt
            fit_results.append({
                'label': f"{cen:.2f} Å",
                'amplitude': amp,
                'center': cen,
                'sigma': sigma,
                'fwhm': 2.355 * sigma,
                'flux': amp * sigma * np.sqrt(2 * np.pi),  # area under Gaussian
                'offset': offset
            })

            # Optional: plot fit
            x_dense = np.linspace(x_fit[0], x_fit[-1], 300)
            plt.plot(x_dense, gaussian(x_dense, *popt), 'r', alpha=0.5)

        except RuntimeError:
            print(f"Could not fit line at ~{x0:.2f} Å")

    sii_region = (real > 6700) & (real < 6745)
    x_fit = real[sii_region]
    y_fit = galnorm[sii_region]

    p0 = [y_fit.max(), 6716.5, 3.0, y_fit.max()/2, np.median(y_fit)]

    try:
        popt, pcov = curve_fit(sii_doublet_with_offset, x_fit, y_fit, p0=p0)

        amp1, cen1, sigma, amp2, offset = popt
        cen2 = cen1 + 14.4

        # Plot the fit
        x_dense = np.linspace(x_fit[0], x_fit[-1], 300)
        plt.plot(x_dense, sii_doublet_with_offset(x_dense, *popt), 'k--')

        print(f"[S II] λ6716: flux = {amp1*sigma*np.sqrt(2*np.pi):.2e}")
        print(f"[S II] λ6731: flux = {amp2*sigma*np.sqrt(2*np.pi):.2e}")

    except RuntimeError:
        print("Fit failed.")

    ha_region = (real > 6525) & (real < 6605)
    x_fit = real[ha_region]
    y_fit = galnorm[ha_region]

    p0 = [y_fit.max(), 6563.0, 3.0, y_fit.max()/5, np.median(y_fit)]

    try:
        popt, pcov = curve_fit(ha_nii_triplet, x_fit, y_fit, p0=p0)
        amp_ha, cen_ha, sigma, amp_nii_6548, offset = popt

        flux_ha = amp_ha * sigma * np.sqrt(2 * np.pi)
        flux_nii_6548 = amp_nii_6548 * sigma * np.sqrt(2 * np.pi)
        flux_nii_6583 = 3 * flux_nii_6548

        print(f"Hα flux: {flux_ha:.2e}")
        print(f"[N II] λ6548 flux: {flux_nii_6548:.2e}")
        print(f"[N II] λ6583 flux: {flux_nii_6583:.2e}")

        x_dense = np.linspace(x_fit[0], x_fit[-1], 300)
        plt.plot(x_dense, ha_nii_triplet(x_dense, *popt), 'k--')

    except RuntimeError:
        print("Hα + [N II] fit failed.")
        
def gaussian(x, amp, cen, sigma, offset):
    return amp * np.exp(-(x - cen)**2 / (2 * sigma**2)) + offset

def sii_doublet(x, amp1, cen1, sigma, amp2):
    cen2 = cen1 + 14.4  # separation in Å between 6716 and 6731
    g1 = amp1 * np.exp(-(x - cen1)**2 / (2 * sigma**2))
    g2 = amp2 * np.exp(-(x - cen2)**2 / (2 * sigma**2))
    return g1 + g2


def sii_doublet_with_offset(x, amp1, cen1, sigma, amp2, offset):
    return sii_doublet(x, amp1, cen1, sigma, amp2) + offset

def ha_nii_triplet(x, amp_ha, cen_ha, sigma, amp_nii_6548, offset):
    cen_nii_6548 = cen_ha - 15.0  # Hα = 6563, [N II] 6548 = -15 Å
    cen_nii_6583 = cen_ha + 20.0  # [N II] 6583 = +20 Å
    amp_nii_6583 = 3 * amp_nii_6548

    g_ha = amp_ha * np.exp(-(x - cen_ha)**2 / (2 * sigma**2))
    g_nii1 = amp_nii_6548 * np.exp(-(x - cen_nii_6548)**2 / (2 * sigma**2))
    g_nii2 = amp_nii_6583 * np.exp(-(x - cen_nii_6583)**2 / (2 * sigma**2))

    return g_ha + g_nii1 + g_nii2 + offset


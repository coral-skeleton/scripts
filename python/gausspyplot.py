import os
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
from pylab import cm
from astropy.wcs import WCS
from astropy import units as u
from gausspyplus.plotting import get_points_for_colormap, shiftedColorMap
import gausspyplus.plotting as gplt
import itertools
import pickle
import random
import sys
import numpy as np
import matplotlib.gridspec as gridspec
from astropy import units as u
from tqdm import tqdm
from gausspyplus.utils.gaussian_functions import gaussian, combined_gaussian
from gausspyplus.utils.spectral_cube_functions import get_spectral_axis, correct_header


#just plot inline
def plot_spectra(pathToDataPickle,
                 path_to_decomp_pickle=None,
                 cols=5,
                 rowsize=7.75,
                 rowbreak=50,
                 dpi=50,
                 n_spectra=None,
                 suffix='',
                 subcube=False,
                 pixel_range=None,
                 list_indices=None,
                 gaussians=True,
                 random_seed=111,
                 vel_unit=u.km/u.s,
                 titlestr = '',
                 label = False,
                 lables = [],
                 combined = False,
                 samecol = False,
                 speccol = False,
                 colors = ['tab:orange'],
                 scale = 1):
    
    
    #get info
    fileName, file_extension = os.path.splitext(os.path.basename(pathToDataPickle))
    data = pickle_load_file(pathToDataPickle)
    decomp = pickle_load_file(path_to_decomp_pickle)
    fileName, file_extension = os.path.splitext(os.path.basename(path_to_decomp_pickle))
    channels, fig_channels = (data['x_values'] for _ in range(2))
    n_channels = len(channels)
    fig_min_channel, fig_max_channel = fig_channels[0], fig_channels[-1] 
    header = gplt.correct_header(data['header'])
    fig_channels = gplt.get_spectral_axis(header=header, to_unit=vel_unit)
    fig_min_channel, fig_max_channel = fig_channels[0], fig_channels[-1]
    list_indices, n_spectra, grid_layout = gplt.get_list_indices(data, subcube=subcube, pixel_range=pixel_range, list_indices=list_indices, n_spectra=n_spectra, random_seed=random_seed)
    cols, rows, rowbreak, colsize, multiple_pdfs = gplt.get_figure_params(n_channels, n_spectra, cols, rowsize, rowbreak, grid_layout, subcube=subcube)
    
    
    #scaling plots to rows
    #fontsize = gplt.scale_fontsize(rowsize)
    figsize = (cols*colsize, rowbreak*rowsize)
    fig = plt.figure(figsize=figsize)

    #actual plot loop
    for i, idx in enumerate(list_indices):

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

        row_i = int((i - (k*rowbreak*cols)) / cols)*3
        col_i = i % cols
        ax = plt.subplot2grid((3*rows_in_figure, cols), (row_i, col_i), rowspan=2)

        ax.step(fig_channels,  y, color = 'black', lw=0.5, where='mid')
        ncomps, rchi2, rchi2gauss, pvalue = None, None, None, None        
        ncomps = len(fit_amps)
        
        if combined:
            combined_gauss = gplt.combined_gaussian(fit_amps, fit_fwhms, fit_means, channels)
            ax.plot(fig_channels, combined_gauss, lw=(3*scale))

            # Plot individual components
        if gaussians:
            for j in range(ncomps):
                gauss = gplt.gaussian(fit_amps[j], fit_fwhms[j], fit_means[j], channels)
                if samecol:
                    ax.plot(fig_channels, gauss, ls='solid', lw=(3*scale), color = colors[0])
                    
                if speccol:
                    if label:
                        ax.plot(fig_channels, gauss, ls='solid', lw=(3*scale), color = colors[j], label = lables[j])
                        ax.legend(fontsize = (14*scale))
                    else:
                        ax.plot(fig_channels, gauss, ls='solid', lw=(3*scale), color = color[j])
                    
                if label:
                    ax.plot(fig_channels, gauss, ls='solid', lw=(3*scale), label = lables[j])
                    ax.legend(fontsize = (14*scale))
                else:
                    ax.plot(fig_channels, gauss, ls='solid', lw=(3*scale))

                
        #axes and title labels            
        xlabel_from_header = r'V$_{los}$ [km/s]'
        btype = 'Intensity '
        bunit = '[Jy/Beam]'
        ylabel_from_header = btype + bunit
        rchi2 = decomp['best_fit_rchi2'][idx]
        loc_string = 'X={}, Y={}'.format(xi, yi)
        title = titlestr
        ax.set_title(title, fontsize=(28*scale))
        ax.set_xlim(fig_min_channel, fig_max_channel)
        ax.set_xlabel(xlabel_from_header, fontsize=(14*scale))
        ax.set_ylabel(ylabel_from_header, fontsize=(14*scale))
        ax.tick_params(labelsize=(12*scale))
        ax.axhline(ls='solid', color = '#191a2b', lw=0.5)

        #row calculations
        if ((i + 1) % (rowbreak*cols) == 0) or ((i + 1) == n_spectra):
            fig.tight_layout()
                    
            remaining_rows = rowbreak

            if (k + 2)*rowbreak > rows:
                remaining_rows = rows - (k + 1)*rowbreak
                figsize = (cols*colsize, remaining_rows*rowsize)

            if (figsize[0] == 0) or (figsize[1] == 0):
                break

            fig = plt.figure(figsize=figsize)
            gridspec.GridSpec(cols, remaining_rows, wspace=0.2, hspace=0.2) 
    
#plot inline with fancy stuff
def plot_spectra_no_save(pathToDataPickle,
                         path_to_decomp_pickle=None,
                         cols=5,
                         rowsize=7.75,
                         rowbreak=50,
                         dpi=50,
                         n_spectra=None,
                         suffix='',
                         subcube=False,
                         pixel_range=None,
                         list_indices=None,
                         gaussians=True,
                         residual=True,
                         signal_ranges=True,
                         random_seed=111,
                         vel_unit=u.km/u.s,
                         label = False, 
                         lables = []):
    
    #load data
    fileName, file_extension = os.path.splitext(os.path.basename(pathToDataPickle))
    data = pickle_load_file(pathToDataPickle)
    decomp = pickle_load_file(path_to_decomp_pickle)
    fileName, file_extension = os.path.splitext(os.path.basename(path_to_decomp_pickle))
    channels, fig_channels = (data['x_values'] for _ in range(2))
    n_channels = len(channels)
    fig_min_channel, fig_max_channel = fig_channels[0], fig_channels[-1] 
    header = gplt.correct_header(data['header'])
    fig_channels = gplt.get_spectral_axis(header=header, to_unit=vel_unit)
    fig_min_channel, fig_max_channel = fig_channels[0], fig_channels[-1]
    list_indices, n_spectra, grid_layout = gplt.get_list_indices(data, subcube=subcube, pixel_range=pixel_range, list_indices=list_indices, n_spectra=n_spectra, 
                                                                 random_seed=random_seed)
    cols, rows, rowbreak, colsize, multiple_pdfs = gplt.get_figure_params(n_channels, n_spectra, cols, rowsize, rowbreak, grid_layout, subcube=subcube)
    
    #scaling
    fontsize = gplt.scale_fontsize(rowsize)
    figsize = (cols*colsize, rowbreak*rowsize)
    fig = plt.figure(figsize=figsize)

    #plot loop
    for i, idx in enumerate(list_indices):
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

        row_i = int((i - (k*rowbreak*cols)) / cols)*3
        col_i = i % cols
        ax = plt.subplot2grid((3*rows_in_figure, cols), (row_i, col_i), rowspan=2)
        if signal_ranges:
            if 'signal_ranges' in data.keys():
                for low, upp in data['signal_ranges'][idx]:
                    ax.axvspan(fig_channels[low], fig_channels[upp - 1], color = '#BBDEF0', alpha = 0.5)
        ax.step(fig_channels,  y, color = 'black', lw=0.5, where='mid')#8aa5b5
        ncomps, rchi2, rchi2gauss, pvalue = None, None, None, None        
        combined_gauss = gplt.combined_gaussian(fit_amps, fit_fwhms, fit_means, channels)
        ncomps = len(fit_amps)
        ax.plot(fig_channels, combined_gauss, lw=2)

        # Plot individual components
        if gaussians:
            for j in range(ncomps):
                gauss = gplt.gaussian(fit_amps[j], fit_fwhms[j], fit_means[j], channels)
                if label:
                    ax.plot(fig_channels, gauss, ls='solid', lw=3, label = lables[j])
                    ax.legend(fontsize = 14)
                else:
                    ax.plot(fig_channels, gauss, ls='solid', lw=3)

                
        ax.plot(fig_channels, combined_gauss, lw=2, color = 'tab:orange')
                    
        
        #axes labels and title            
        xlabel_from_header = r'V$_{los}$ [km/s]'
        if 'CTYPE3' in header.keys():
            xlabel_from_header = '{} [{}]'.format(header['CTYPE3'], vel_unit)
        btype = 'Intensity '
        bunit = '[Jy/Beam]'
        ylabel_from_header = btype + bunit
        rchi2 = decomp['best_fit_rchi2'][idx]
        idx_string = ' (Idx$_{{data}}$={})'.format(index_data)
        loc_string = 'X={}, Y={}'.format(xi, yi)
        ncomps_string = ', N$_{{comp}}$={}'.format(ncomps)
        rchi2_string = ', $\\chi_{{red}}^{{2}}$={:.3f}'.format(rchi2)
        title = 'Idx={}{}{}{}{}'.format(idx, idx_string, loc_string, ncomps_string, rchi2_string)
        title = '{}'.format(loc_string)
        ax.set_title(title, fontsize=28)
        ax.set_xlim(fig_min_channel, fig_max_channel)
        ax.set_xlabel(xlabel_from_header, fontsize=14)
        ax.set_ylabel(ylabel_from_header, fontsize=14)
        ax.tick_params(labelsize=12)
        ax.axhline(ls='solid', color = '#191a2b', lw=0.5)
        ax.axhline(y=rms, color = '#303354', ls='dotted', lw=0.5)
        ax.axhline(y=-rms, color = '#303354',  ls='dotted', lw=0.5)
        ax.axhline(y=3*rms, color = '#f71b1b', ls='dashed', lw=1)
    
        #row loop
        if ((i + 1) % (rowbreak*cols) == 0) or ((i + 1) == n_spectra):
            fig.tight_layout()
                    
            remaining_rows = rowbreak

            if (k + 2)*rowbreak > rows:
                remaining_rows = rows - (k + 1)*rowbreak
                figsize = (cols*colsize, remaining_rows*rowsize)

            if (figsize[0] == 0) or (figsize[1] == 0):
                break

            fig = plt.figure(figsize=figsize)
            gridspec.GridSpec(cols, remaining_rows, wspace=0.2, hspace=0.2)     
            
def get_cmap_rchi2(vmin, vmax):
    orig_cmap = matplotlib.cm.RdBu_r
    start, stop = get_points_for_colormap(vmin, vmax, central_val=1.)
    midpoint = (1 - vmin) / (vmax - vmin)
    return shiftedColorMap(orig_cmap, start=0., midpoint=midpoint, stop=stop)

def add_style(ax):
    ax.set_xlabel('Galactic Longitude')
    ax.set_ylabel('Galactic Latitude')

def get_points_for_colormap(vmin, vmax, central_val=0.):
    lower_interval = abs(central_val - vmin)
    upper_interval = abs(vmax - central_val)

    if lower_interval > upper_interval:
        start = 0.
        stop = 0.5 + (upper_interval / lower_interval)*0.5
    else:
        start = 0.5 - (lower_interval / upper_interval)*0.5
        stop = 1.
    return start, stop

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def pickle_load_file(pathToFile):
    with open(os.path.join(pathToFile), "rb") as pickle_file:
        if (sys.version_info > (3, 0)):
            data = pickle.load(pickle_file, encoding='latin1')
        else:
            data = pickle.load(pickle_file)
    return data

def get_list_indices(data, subcube=False, pixel_range=None,
                     list_indices=None, n_spectra=None, random_seed=111):
    random.seed(random_seed)
    # TODO: incorporate the nan_mask in this scheme
    grid_layout = None
    if subcube or (pixel_range is not None):
        if pixel_range is not None:
            xmin, xmax = pixel_range['x']
            ymin, ymax = pixel_range['y']
        else:
            ymin, xmin = min(data['location'])
            ymax, xmax = max(data['location'])
        yValues = np.arange(ymin, ymax + 1)[::-1]
        xValues = np.arange(xmin, xmax + 1)
        locations = list(itertools.product(yValues, xValues))
        cols = len(xValues)
        rows = len(yValues)
        grid_layout = [cols, rows]
        n_spectra = cols*rows

        list_indices = []
        for location in locations:
            list_indices.append(data['location'].index(location))
    elif list_indices is None:
        if n_spectra is None:
            n_spectra = len(data['data_list'])
            list_indices = np.arange(n_spectra)
        else:
            list_indices = []
            nIndices = len(data['data_list'])
            randomIndices = random.sample(range(nIndices), nIndices)
            for idx in randomIndices:
                if 'nan_mask' in data.keys():
                    yi, xi = data['location'][idx]
                    if data['nan_mask'][:, yi, xi].all() != True:
                        list_indices.append(idx)
                        if len(list_indices) == n_spectra:
                            break
                else:
                    list_indices.append(idx)
                    if len(list_indices) == n_spectra:
                        break
    else:
        n_spectra = len(list_indices)

    list_indices = [i for i in list_indices if data['data_list'][i] is not None]
    n_spectra = len(list_indices)

    return list_indices, n_spectra, grid_layout

def get_figure_params(n_channels, n_spectra, cols, rowsize, rowbreak,
                      grid_layout, subcube=False):
    if n_channels > 700:
        colsize = round(rowsize*n_channels/659, 2)
    else:
        colsize = rowsize

    # if subcube is False:
    if grid_layout is None:
        rows = int(n_spectra / (cols))
        if n_spectra % cols != 0:
            rows += 1

        multiple_pdfs = True
        if rows < rowbreak:
            rowbreak = rows
            multiple_pdfs = False
    else:
        cols, rows = grid_layout
        rowbreak = rows
        multiple_pdfs = False

    if (rowbreak*rowsize*100 > 2**16) or (cols*colsize*100 > 2**16):
        errorMessage = \
            "Image size is too large. It must be less than 2^16 pixels in each direction. Restrict the number of columns or rows."
        raise Exception(errorMessage)

    return cols, rows, rowbreak, colsize, multiple_pdfs

def xlabel_from_header(header, vel_unit):
    xlabel = 'Channels'

    if header is None:
        return xlabel

    if 'CTYPE3' in header.keys():
        xlabel = '{} [{}]'.format(header['CTYPE3'], vel_unit)

    return xlabel

def ylabel_from_header(header):
    if header is None:
        return 'Intensity'

    btype = 'Intensity'
    if 'BTYPE' in header.keys():
        btype = header['BTYPE']

    bunit = ''
    if 'BUNIT' in header.keys():
        bunit = ' [{}]'.format(header['BUNIT'])
        
    bunit = ' [Jy/Beam]'
    return btype + bunit

def add_figure_properties(ax, rms, fig_min_channel, fig_max_channel,
                          header=None, residual=False, fontsize=10, vel_unit=u.km/u.s):
    # TODO: read labels automatically from header
    ax.set_xlim(fig_min_channel, fig_max_channel)
    ax.set_xlabel(xlabel_from_header(header, vel_unit), fontsize=fontsize)
    ax.set_ylabel(ylabel_from_header(header), fontsize=fontsize)

    ax.tick_params(labelsize=fontsize - 2)

    ax.axhline(color='black', ls='solid', lw=0.5)
    ax.axhline(y=rms, color='red', ls='dotted', lw=0.5)
    ax.axhline(y=-rms, color='red', ls='dotted', lw=0.5)

    if not residual:
        ax.axhline(y=3*rms, color='red', ls='dashed', lw=1)
    else:
        ax.set_title('Residual', fontsize=fontsize)
        
def plot_signal_ranges(ax, data, idx, fig_channels):
    if 'signal_ranges' in data.keys():
        for low, upp in data['signal_ranges'][idx]:
            ax.axvspan(fig_channels[low], fig_channels[upp - 1], alpha=0.1, color='indianred')

            
def get_title_string(idx, index_data, xi, yi, ncomps, rchi2, rchi2gauss, pvalue):
    idx_string = ''
    if index_data is not None:
        if index_data != idx:
            idx_string = ' (Idx$_{{data}}$={})'.format(index_data)

    loc_string = ''
    if xi is not None:
        loc_string = ', X={}, Y={}'.format(xi, yi)

    ncomps_string = ''
    if ncomps is not None:
        ncomps_string = ', N$_{{comp}}$={}'.format(ncomps)

    rchi2_string = ''
    if rchi2 is not None:
        rchi2_string = ', $\\chi_{{red}}^{{2}}$={:.3f}'.format(rchi2)

    rchi2gauss_string = ''
    if rchi2gauss is not None:
        rchi2gauss_string = ', $\\chi_{{red, gauss}}^{{2}}$={:.3f}'.format(
            rchi2gauss)

    pvalue_string = ''
    # if pvalue is not None:
    #     pvalue_string = ', $p_{{value}}$={:.3f}'.format(
    #         pvalue)

    title = 'Idx={}{}{}{}{}{}{}'.format(
        idx, idx_string, loc_string, ncomps_string, rchi2_string,
        rchi2gauss_string, pvalue_string)
    return title

def scale_fontsize(rowsize):
    rowsize_scale = 4
    if rowsize >= rowsize_scale:
        fontsize = 10 + int(rowsize - rowsize_scale)
    else:
        fontsize = 10 - int(rowsize - rowsize_scale)
    return fontsize

def plot_spectra_save(pathToDataPickle, *args,
                 path_to_plots=None, path_to_decomp_pickle=None,
                 training_set=False, cols=5, rowsize=7.75, rowbreak=50, dpi=5000,
                 n_spectra=None, suffix='', subcube=False, pixel_range=None,
                 list_indices=None, gaussians=True, residual=True,
                 signal_ranges=True, random_seed=111, vel_unit=u.km/u.s):

    print("\nPlotting...")

    #  for compatibility with older versions
    if args:
        path_to_plots = args[0]

    if path_to_plots is None:
        if path_to_decomp_pickle is None:
            path_to_plots = os.path.dirname(pathToDataPickle)
        else:
            path_to_plots = os.path.dirname(path_to_decomp_pickle)

    decomposition = True
    if path_to_decomp_pickle is None:
        decomposition = False

    fileName, file_extension = os.path.splitext(os.path.basename(pathToDataPickle))

    data = pickle_load_file(pathToDataPickle)
    if (not training_set) and decomposition:
        decomp = pickle_load_file(path_to_decomp_pickle)
        fileName, file_extension = os.path.splitext(
            os.path.basename(path_to_decomp_pickle))

    channels, fig_channels = (data['x_values'] for _ in range(2))
    n_channels = len(channels)
    fig_min_channel, fig_max_channel = fig_channels[0], fig_channels[-1]

    if 'header' in data.keys():
        header = correct_header(data['header'])
        fig_channels = get_spectral_axis(header=header, to_unit=vel_unit)
        fig_min_channel, fig_max_channel = fig_channels[0], fig_channels[-1]
    else:
        header = None

    list_indices, n_spectra, grid_layout = get_list_indices(
        data, subcube=subcube, pixel_range=pixel_range, list_indices=list_indices, n_spectra=n_spectra, random_seed=random_seed)

    cols, rows, rowbreak, colsize, multiple_pdfs = get_figure_params(
        n_channels, n_spectra, cols, rowsize, rowbreak, grid_layout, subcube=subcube)

    fontsize = scale_fontsize(rowsize)

    figsize = (cols*colsize, rowbreak*rowsize)
    fig = plt.figure(figsize=figsize)

    fig.patch.set_facecolor('white')
    fig.patch.set_alpha(1.0)
    # fig.subplots_adjust(hspace=0.5)

    pbar = tqdm(total=n_spectra)

    for i, idx in enumerate(list_indices):
        pbar.update(1)
        if 'location' in data.keys():
            yi, xi = data['location'][idx]
        else:
            yi, xi = (None for _ in range(2))

        if 'index' in data.keys():
            index_data = data['index'][idx]
        else:
            index_data = None

        if rowbreak is not None:
            k = int(i / (rowbreak*cols))
            if (k + 1)*rowbreak > rows:
                rows_in_figure = rows - k*rowbreak
            else:
                rows_in_figure = rowbreak

        y = data['data_list'][idx]
        rms = data['error'][idx][0]

        if training_set:
            fit_fwhms = data['fwhms'][idx]
            fit_means = data['means'][idx]
            fit_amps = data['amplitudes'][idx]
        elif decomposition:
            fit_fwhms = decomp['fwhms_fit'][idx]
            fit_means = decomp['means_fit'][idx]
            fit_amps = decomp['amplitudes_fit'][idx]

        row_i = int((i - (k*rowbreak*cols)) / cols)*3
        col_i = i % cols
        ax = plt.subplot2grid((3*rows_in_figure, cols),
                              (row_i, col_i), rowspan=2)

        ax.step(fig_channels, y, color='black', lw=0.5, where='mid')

        ncomps, rchi2, rchi2gauss, pvalue = None, None, None, None

        if decomposition or training_set:
            combined_gauss = combined_gaussian(
                fit_amps, fit_fwhms, fit_means, channels)
            ax.plot(fig_channels, combined_gauss, lw=2, color='orangered')

            ncomps = len(fit_amps)

            # Plot individual components
            if gaussians:
                for j in range(ncomps):
                    gauss = gaussian(
                        fit_amps[j], fit_fwhms[j], fit_means[j], channels)
                    ax.plot(fig_channels, gauss, ls='solid', lw=1, color='orangered')

        if training_set:
            rchi2 = data['best_fit_rchi2'][idx]
            if 'pvalue' in data.keys():
                pvalue = data['pvalue'][idx]

        elif decomposition:
            rchi2 = decomp['best_fit_rchi2'][idx]

        if signal_ranges:
            plot_signal_ranges(ax, data, idx, fig_channels)

        rchi2gauss = None
        # TODO: incorporate rchi2_gauss

        title = get_title_string(idx, index_data, xi, yi, ncomps,
                                 rchi2, rchi2gauss, pvalue)
        ax.set_title(title, fontsize=fontsize)

        add_figure_properties(ax, rms, fig_min_channel, fig_max_channel,
                              header=header, fontsize=fontsize, vel_unit=vel_unit)

        if residual and (decomposition or training_set):
            row_i = int((i - k*(rowbreak*cols)) / cols)*3 + 2
            col_i = i % cols
            ax = plt.subplot2grid((3*rows_in_figure, cols),
                                  (row_i, col_i))

            ax.step(fig_channels, y - combined_gauss, color='black', lw=0.5, where='mid')
            if signal_ranges:
                plot_signal_ranges(ax, data, idx, fig_channels)

            add_figure_properties(ax, rms, fig_min_channel, fig_max_channel, header=header,
                                  residual=True, fontsize=fontsize, vel_unit=vel_unit)

        if ((i + 1) % (rowbreak*cols) == 0) or ((i + 1) == n_spectra):
            if multiple_pdfs:
                filename = '{}{}_plots_part_{}.pdf'.format(fileName, suffix, k + 1)
            else:
                filename = '{}{}_plots.pdf'.format(fileName, suffix)

            fig.tight_layout()

            if not os.path.exists(path_to_plots):
                os.makedirs(path_to_plots)
            pathname = os.path.join(path_to_plots, filename)
            fig.savefig(pathname, dpi=dpi, overwrite=True)
            plt.close()

            #  close progress bar before print statement to avoid duplicate
            #  progress bars
            if pbar.n >= n_spectra:
                pbar.close()
            print("\n\033[92mSAVED FILE:\033[0m '{}' in '{}'".format(filename, path_to_plots))

            remaining_rows = rowbreak
            #  for last iteration
            if (k + 2)*rowbreak > rows:
                remaining_rows = rows - (k + 1)*rowbreak
                figsize = (cols*colsize, remaining_rows*rowsize)

            if (figsize[0] == 0) or (figsize[1] == 0):
                break

            fig = plt.figure(figsize=figsize)

            # set up new subplot grid
            gridspec.GridSpec(cols, remaining_rows, wspace=0.2, hspace=0.2)

            fig.patch.set_facecolor('white')
            fig.patch.set_alpha(1.0)
            fig.subplots_adjust(hspace=0.5)
    plt.close()
    

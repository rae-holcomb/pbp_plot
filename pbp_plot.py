import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
from scipy import stats

# number of seconds in a day
day = 24*60*60

def download_tpf(target_name, sector=None) :
    """
    Given the target name, i.e. 'TIC 233629963', 'TOI 1468', or 'LSPM J1835+3259', downloads the target pixel
    file for the first available sector. If sector keyword is specified (as an integer), that sector only will
    be retrieved. Multisector downloads are not currently supported.
        Input:  target_name - string
                sector - integer of which TESS sector you want
        Output: tpf - a single sector lightkurve TargetPixelFile object 
    """
    if sector :
        return lk.search_targetpixelfile(target_name, sector=sector).download()
    else :
        return lk.search_targetpixelfile(target_name).download()


def plot_pixel_by_pixel(target_name, tpf, plot_type='lc', padding=1, x_lim=None, y_lim=None, bin=False, save=False, figname=None) :
    """
    Creates a pixel-by-pixel plot of a tess target from the target pixel file. Subplots
    in blue are pixels inside the TESS pipeline aperture and those in red are outside it.
    Note that this code is best used for investigating signals of know duration or period,
    the plots are small so it is very helpful to know your x_lim before plotting. Note that
    this may run slowly, especially if a large aperture or padding is used.
    
    Inputs: target_name - the name of the target to be plotted. Must be resolvable by the TIC.
            plot_type - either 'lc' or 'periodogram'. You may also create a custom function with a bit more work.
            padding - size of the pixel border around the TESS aperture you want to plot. Note
                      that a large padding can result in slow runtime to generate plots.
            x_lim - for 'lc' the BJD dates you want plotted. For 'periodogram' the range of period in days.
            y_lim - the limits of the y-axis. If left blank with be scaled dynamically.
            save - Boolean of whether you would like the plots saved to a file or not.
    """
    def plot_lc(tpf, mask, ax=None, c='b', x_lim=x_lim, y_lim=y_lim, bin=bin) :
        """
        Given a masked target pixel file, calculates the light curve and returns it to be
        plotted in a subplot of the pixel-by-pixel plot.
            Input:  tpf - A target pixel file
                    mask - the associated aperture mask for the tpf
                    ax - the axes object to be plotted into
                    c - the color of the plot
                    bin - if set to a float, bins the lc to time bins of that size in units of days. may run slowly
            Output: x, y - arrays to be plotted
        """
        ax = ax or plt.gca()

        # calculate light curve
        lc = tpf.to_lightcurve(aperture_mask=mask.astype(bool)).flatten()
        # bin_lc = lc.flatten(window_length=801).remove_outliers().bin(binsize=5)

        # bin if needed
        if bin:
            time, flux = bin_lc(lc.time.value, lc.flux.value, bs=bin*60*60*24)
        else:
            time = lc.time.value
            flux = lc.flux.value

        # do the plotting
        ax.plot(time,flux, linestyle='', marker='.', c=c)

        # set x axis window (currently in days)
        if x_lim :
            ax.set_xlim(x_lim)
        else :
            x_lim = [lc.time.value[0], lc.time.value[-1]]
            ax.set_xlim(x_lim)
        
        # set y axis window
        x_mask = np.where((lc.time.value > x_lim[0]) & (lc.time.value < x_lim[1]))
        if y_lim :
            ax.set_ylim(y_lim)
        else :
            ax.set_ylim([np.nanpercentile(lc.flux[x_mask], 5) * 0.8, np.nanpercentile(lc.flux[x_mask], 97) * 1.2])

        return ax

    def plot_periodogram(tpf, mask, ax=None, c='b', x_lim=x_lim, y_lim=y_lim, bin=bin):
        """
        Given a masked target pixel file, calculates the periodogram and returns it to be
        plotted in a subplot of the pixel-by-pixel plot.
            Input:  tpf - A target pixel file
                    mask - the associated aperture mask for the tpf
                    ax - the axes object to be plotted into
                    c - the color of the plot
                    Output: x, y - arrays to be plotted
        """    
        ax = ax or plt.gca()

        # calculate periodogram
        lc = tpf.to_lightcurve(aperture_mask=mask.astype(bool)).flatten()
        periodogram = lc.to_periodogram()
        periodogram.default_view='period'

        # plot
        line = ax.plot(periodogram.period, periodogram.power, c=c)

        # set x axis window (currently in days)
        if x_lim :
            ax.set_xlim(x_lim)
        else :
            x_lim = [.1,10]
            ax.set_xlim(x_lim)

        # fix the scaling of the y-axis
        if y_lim :
            ax.set_ylim(y_lim)
        else : 
            y_window = np.where((np.array(periodogram.period) < x_lim[1]) & (np.array(periodogram.period) > x_lim[0]))[0]
            max_power = np.max(np.array(periodogram.power)[y_window])
            y_lim = [0, max_power*1.2]
            ax.set_ylim(y_lim)

        return line

    # feel free to add in whatever plotting function you like here
    def plot_custom(tpf, mask, ax=None, c='b', x_lim=x_lim, y_lim=y_lim) :
        """
        To be populated by user.
        Given a masked target pixel file, calculates a custom plot and returns it to be
        plotted in a subplot of the pixel-by-pixel plot.
            Input:  tpf - A target pixel file
                    mask - the associated aperture mask for the tpf
                    Output: x, y - arrays to be plotted
        """
        ax = ax or plt.gca()
        # do what you will with the tpf
        return

    
    # Make plot of the TESS pipeline aperture
    tess_fig = tpf.plot(aperture_mask=tpf.pipeline_mask)
    if save :
        plot_name = target_name + '_tpf_figure'
        plt.savefig(plot_name, facecolor=tess_fig.get_facecolor(), edgecolor='none')
    else :
        plt.show()
        
    # Check that the provided x_lim values are in range. If not, fix them.
    if plot_type == 'lc' :
        if not x_lim :
            x_lim = [tpf[0].time.value, tpf[-1].time.value]
        if x_lim[0] < tpf[0].time.value :
            print('x_lim should be provided in BJD - 2457000 (BTJD). The valid date range for this system is ' + str([tpf[0].time.value, tpf[-1].time.value]) + '.')
            x_lim[0] = tpf[0].time.value
        if x_lim[1] > tpf[-1].time.value :
            print('x_lim should be provided in BJD - 2457000 (BTJD). The valid date range for this system is ' + str([tpf[0].time.value, tpf[-1].time.value]) + '.')
            x_lim[1] = tpf[-1].time.value
    elif plot_type == 'periodogram' :
        if not x_lim :
            print('x_lim for periodograms should be provided as a [min_period, max_period] in units of days.')
            x_lim = [.1, 8]
        elif x_lim[0] < 0 :
            print('x_lim for periodograms should be provided as a [min_period, max_period] in units of days.')
            xlim[0] = 0

    
    # Create our custom mask
    pipe_mask = tpf.pipeline_mask

    # This grabs the dimensions of the existing pipeline mask
    startrow = np.where(np.sum(pipe_mask, axis=1) > 0)[0][0]
    endrow = np.where(np.sum(pipe_mask, axis=1) > 0)[0][-1] + 1
    startcol = np.where(np.sum(pipe_mask, axis=0) > 0)[0][0]
    endcol = np.where(np.sum(pipe_mask, axis=0) > 0)[0][-1] + 1

    # Get the dimensions of an aperature of the size of the pipeline mask 
    # with extra padding around it
    # (abbr for e.g. "start row expanded")
    if (startrow - padding > 0): sre = startrow - padding
    else: sre = 0
    if (startcol - padding > 0): sce = startcol - padding
    else: sce = 0
    if (endrow + padding < tpf.shape[1:][0]): ere = endrow + padding
    else: ere = tpf.shape[1:][0]
    if (endcol + padding < tpf.shape[1:][1]): ece = endcol + padding
    else: ece = tpf.shape[1:][1]

    # create custom mask
    custom_apt = np.zeros(tpf.shape[1:], dtype=np.int)
    custom_apt[sre:ere, sce:ece] = 1

    # set up the subplots
    fig, axs = plt.subplots(ere-sre, ece-sce, figsize=[18,15])
    fig.set_facecolor('white')
    plt.tight_layout()

    # add overall axes labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.grid(False)
    if figname :
        plt.title(figname, fontsize=26)
    else: 
        plt.title(target_name, fontsize=26)

    # set the plot function and label axes correctly.
    if plot_type == 'lc' :
        plot_function = plot_lc
        plt.xlabel("Time (BJD)",fontsize=18)
        plt.ylabel("Normalized Flux",fontsize=18)

    elif plot_type == 'periodogram' :
        plot_function = plot_periodogram
        plt.xlabel("Period (days)",fontsize=18)
        plt.ylabel("LSP Power",fontsize=18)

    else :
        print("Invalid plot type. Try 'lc' or 'periodogram'.")
        

    # loop through each pixel in frame
    for i in range(sre, ere) :
        for j in range(sce,ece) : 
            # make aperture
            apij = np.zeros(tpf.shape[1:], dtype=int)
            apij[i][j] = 1

            # plot the line
            # check if it is included in the mask and color accordingly
            if tpf.pipeline_mask[i][j] :
                ax = plot_function(tpf, apij, ax=axs[(ere-i-1), j-sce], c='b')
            else :
                ax = plot_function(tpf, apij, ax=axs[(ere-i-1), j-sce], c='r')

    plt.tight_layout()
    
    if save :
        plot_name = target_name + '_PbP_' + plot_type
        fig.savefig(plot_name, facecolor='white', edgecolor='none')
    else :
        plt.show()
 
def bin_lc(time, flux, flux_err=None, bs=120):
    """
    Bins a timeseries to the desired cadence. Works much faster than Lightkurve's built in binning function.

    Args:
        time (:obj:`array`): time stamps of the lightcurve.
        flux (:obj:`array`): flux values of the lightcurve.
        flux_err (:obj:`array`): flux error values of the lightcurve.
        bs (:obj:`float`): the size of the bins, in units of seconds.

    Returns:
    Two parameters, or three if a flux_err is also provided.

        - time_bin (:obj:`array`): binned time stamps of the lightcurve.
        - flux_bin (:obj:`array`): binned flux values of the lightcurve.
        - flux_err_bin (:obj:`array`): binned flux errors values of the lightcurve. Only returned if an array is passed to flux_err.
    """
    time_bin = np.arange(time[0], time[-1], bs/day)
    flux_bin = stats.binned_statistic(time, flux, bins=time_bin)[0]

    # If flux_err is populated, assume the errors combine as the root-mean-square
    if flux_err is not None:
        # define a function to calculate the root mean square error of each bin
        rmse_func = (lambda x: np.sqrt(np.nansum(np.square(x))) / len(np.atleast_1d(x))
                    if np.any(np.isfinite(x))
                    else np.nan)
        flux_err_bin = stats.binned_statistic(time, flux_err, statistic=rmse_func, bins=time_bin)[0]
        return time_bin[:-1], flux_bin, flux_err_bin

    return time_bin[:-1], flux_bin
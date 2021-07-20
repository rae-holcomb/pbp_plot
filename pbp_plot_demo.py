import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk
import pbp_plot


# The name target you would like to plot
target_name = 'LSPM J1835+3259'

# Download the target pixel file
tpf = pbp_plot.download_tpf(target_name)

# Set the paramters to plot a light curve
x_lim = [2016.75,2017.5]	# time range in BJD - 2457000 (BTJD)
padding = 1
pbp_plot.plot_pixel_by_pixel(target_name, tpf, plot_type='lc', padding=1, x_lim=x_lim, save=True)


# Set the paramters to plot a periodogram
period_lim = [.5,4.]	# range of periods in units of days
padding = 1
pbp_plot.plot_pixel_by_pixel(target_name, tpf, plot_type='periodogram', padding=1, x_lim=period_lim, save=True)

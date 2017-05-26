import numpy as np
from waveletFunctions import wavelet, wave_signif
import matplotlib.pylab as plt
import matplotlib.dates as dates
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from sunpy import lightcurve as lc
from sunpy.time import parse_time
from sunpy import lightcurve as lc
import datetime
import scipy
from matplotlib import mlab, cm
from pandas import*
from sunpy.instr import lyra
from astropy.io import fits
import datetime
from scipy.io.idl import readsav
smooth = scipy.ndimage.filters.uniform_filter
plt.ion()
from scipy import fftpack
from sunpy.time import parse_time
basetime=parse_time('1979-01-01 00:00:00')
import datetime


import seaborn as sns
sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
sns.set_context('paper', font_scale = 1.5)

'''
    
Function to plot wavelet analysis.
Takes np array, and dt.
'''

def wave_fn(sst, dt, title = 'test', cmapp = 'Greys'):
    
    
    variance = np.std(sst, ddof=1) ** 2
    sst = (sst - np.mean(sst)) / np.std(sst, ddof=1)
    n = len(sst)
    
    xlim = ([0, dt*n])
    time = np.arange(len(sst)) * dt  # construct time array
    
    pad = 1  # pad the time series with zeroes (recommended)
    dj = 0.25  # this will do 4 sub-octaves per octave
    s0 = 2 * dt  # this says start at a scale of 6 months
    j1 = 7 / dj  # this says do 7 powers-of-two with dj sub-octaves each
    #lag1 = 0.72  # lag-1 autocorrelation for red noise background
    lag1 = ((np.corrcoef(sst, np.roll(sst,1)) + np.sqrt(np.corrcoef(sst, np.roll(sst, 2))))/2)[1][0]
    mother = 'MORLET'
    
    # Wavelet transform:
    wave, period, scale, coi = wavelet(sst, dt, pad, dj, s0, j1, mother)
    powers = (np.abs(wave)) ** 2  # compute wavelet power spectrum
    power=np.zeros_like(powers)
    for k in range(len(scale)):
        power[k,:] = powers[k,:]/period[k]


    # Significance levels: (variance=1 for the normalized SST)
    signif = wave_signif(([1.0]), dt=dt, sigtest=0, scale=scale, lag1 = lag1, mother=mother)
    sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand signif --> (J+1)x(N) array
    sig95 = power / sig95  # where ratio > 1, power is significant
    
    # Global wavelet spectrum & significance levels:
    global_ws = variance * (np.sum(power, axis=1) / n)  # time-average over all times
    dof = n - scale  # the -scale corrects for padding at edges
    global_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=1, lag1 = lag1, dof=dof, mother=mother)
    
    # Scale-average between El Nino periods of 2--8 years
    avg = np.logical_and(scale >= 2, scale < 8)
    Cdelta = 0.776  # this is for the MORLET wavelet
    scale_avg = scale[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])  # expand scale --> (J+1)x(N) array
    scale_avg = power / scale_avg  # [Eqn(24)]
    scale_avg = variance * dj * dt / Cdelta * sum(scale_avg[avg, :])  # [Eqn(24)]
    scaleavg_signif = wave_signif(variance, dt=dt, scale=scale, sigtest=2, dof=([2, 7.9]), mother=mother)
    
    
    ##################################################
    #-------------------Plotting---------------------#
    ##################################################
    
    
    #--- Plot time series
    plt.figure(figsize=(7, 8))
    plt.subplot(211)
    
    #plt.plot(time,smooth(gl_deriv,10))
    #plt.plot(time, gl_1)
    #plt.plot(time, gl_2)
    plt.plot(time, sst)
    plt.xlim(xlim[:])
    plt.ylim(np.min(sst), np.max(sst))
    #plt.ylim(np.min(smooth(gl_deriv,10)), np.max(smooth(gl_deriv,5)))
    plt.tick_params(
                    axis='x',          # changes apply to the x-axis
                    which='both',      # both major and minor ticks are affected
                    bottom='off',      # ticks along the bottom edge are off
                    top='off',         # ticks along the top edge are off
                    labelbottom='off')
                    #plt.xlabel('Time (seconds)')
    plt.ylabel('Detrended Lighcurve')
    plt.title(title)
    plt.hold(False)
    
    
    #--- Contour plot wavelet power spectrum
    plt3 = plt.subplot(212)
    
    
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16]
    #levels = [0, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.6]
    norm = cm.colors.Normalize(vmax=abs(power).max(), vmin=np.min(power))
    levels = np.linspace(0, np.max(power), 200)
    cmap = cmapp
    
    CS = plt.contour(time, scale, (power), len(levels))  #*** or use 'contour'
    im = plt.contourf(CS, levels = levels, cmap = cmap, norm = norm)#, norm = LogNorm())
    plt.xlabel('X values')
    plt.ylabel('Wavelet Time Scale')
    
    plt3.invert_yaxis()
    #plt.title('b) NINO3 SST Wavelet Power Spectrum (in base 2 logarithm)')
    plt.xlim(xlim[:])
    
        
        
        

    # 95# significance contour, levels at -99 (fake) and 1 (95# signif)
    plt.hold(True)
    #plt.contour(time, scale, sig95, [-99, 1], colors='white')
    # cone-of-influence, anything "below" is dubious
    plt.plot(time, coi, 'white')
    #plt.hold(False)
    # format y-scale
    #plt3.set_yscale('log', basey=2, subsy=None)
    plt.ylim(np.min(period), period[len(period)-2])
    ax = plt.gca().yaxis
    ax.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    plt3.ticklabel_format(axis='y', style='plain')
    #plt3.invert_yaxis()
    # set up the size and location of the colorbar
    divider = make_axes_locatable(plt3)
    #cax = divider.append_axes("bottom", size="5%", pad=0.5)
    #plt.colorbar(im, cax=cax, orientation='horizontal')
    
    plt.tight_layout()
    plt.subplots_adjust(hspace = 0.03)
    plt.show()

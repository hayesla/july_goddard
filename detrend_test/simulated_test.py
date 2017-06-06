import numpy as np
import matplotlib.pyplot as plt
import scipy
from pandas import*
import datetime
from scipy import ndimage
smooth = ndimage.filters.uniform_filter
plt.ion()
from scipy import fftpack
from scipy.signal import savgol_filter
import os
import astropy.units as u
from timeseries_simulation import TimeSeriesFromModelSpectrum
from rnspectralmodels4 import power_law
from scipy.special import erf
from sunpy import lightcurve as lc
from simulated_flare import SimulatedFlare
from scipy.signal import chirp


import seaborn as sns
sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
sns.set_context('paper', font_scale = 1.5)


def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def normalise(x): #Function to normalise data
	return (x-x.min())/(x.max() - x.min())

def phi(x):
    return (1./(np.sqrt(2*np.pi)))*np.exp(-np.power(x, 2)/2)

def PHI(x):
    return (1./2)*(1 + erf(x/np.sqrt(2)))


t_start = '2012-07-19 04:15:00'
t_end = '2012-07-19 09:27:00'

#flare background
goes_lc = lc.GOESLightCurve.create(t_start, t_end)
goeslong = goes_lc.data['xrsb']
test_lc = normalise(np.array(goeslong))
x = np.arange(0, len(test_lc))
z = np.polyfit(x, test_lc, 20)
p = np.poly1d(z)
flare_background = smooth(p(x),1000)
n = len(flare_background)
dt = 1.

#white noise
white_noise = np.random.normal(0, 0.0001, n)

#red noise
r_n = TimeSeriesFromModelSpectrum(power_law, [1, 2], n, dt, 1, 1)
red_noise = np.real(r_n.sample())/10

#increaing period sin wave
period = np.linspace(10, 50, len(x))
freq = 1./period
#freq = 1./10
amp = np.linspace(0.001, 0.0001, len(x))
yy = test_lc - smooth(test_lc, 70)
amp = smooth(np.abs(yy), 70)
#amp = np.logspace(0.001, 0.0005, len(x))
#amp = gaussian(x, 2000, 800) + 1
sin_com = []
for i in range(len(x)):
    sin_com.append(amp[i]*np.sin(2*np.pi*freq[i]*x[i]))
sin_com = np.array(sin_com)


background1 = np.zeros(len(x))
ts = SimulatedFlare(flare_background, white_noise, red_noise, sin_com, n, dt)

ts1 = SimulatedFlare(background1, white_noise, red_noise, sin_com, n, dt)

from wavelet_test import wave_fn


def fourier_ana(x):
    N = len(x)
    dt = 2.0
    df = 1./(N*dt)
    PSD = abs(dt*fftpack.fft(x)[:N/2])**2
    f = df*np.arange(N/2)
    fp = f > 1./1000 #only want periods < 1000
    f = f[fp]
    PSD = PSD[fp]
    #plt.plot(1./f, PSD/np.max(PSD))
    return f, PSD



save_dir = '/Users/laura/Documents/QPPs/July_event/july_goddard/detrend_test/smooth_lcs'
def simulated_data_check():
    smooth_factor = np.arange(11, 301, 6)
    smooth_factor_sg = np.arange(11, 301, 6)
    a, b = fourier_ana(sin_com)
    data_cube = []
    data_cube2 = []
    for i in range(len(smooth_factor)):
        data_cube.append(ts.flare - savgol_filter(ts.flare, smooth_factor[i], 3))
        data_cube2.append(ts.flare - smooth(ts.flare, smooth_factor[i]))
    for i in range(len(data_cube)):
        f, PSD = fourier_ana(data_cube[i])
        f1, PSD1 = fourier_ana(data_cube2[i])
        fig, ax = plt.subplots(3, figsize = (7,8) )
        ax[0].plot(data_cube[i]/np.max(data_cube[i]), label = 'Savgol window = ' + str(smooth_factor[i]) + ' s', color = 'g')
        
        ax[0].grid()
        ax[0].legend()
        ax[0].set_ylabel('Detrended Flux')
        
        ax[1].plot(data_cube2[i]/np.max(data_cube2[i]), label = 'Smooth window = '+str(smooth_factor[i]), color = 'r')
        ax[1].grid()
        ax[1].legend()
        ax[1].set_ylabel('Detrended Flux')
        
        ax[2].plot(1./f, PSD,label = 'FT', color = 'g')
        ax[2].plot(1./f1, PSD1,label = 'FT', color = 'r')
        ax[2].plot(1./a, b, label = 'Actual Period', color = 'k')
        ax[2].set_xlabel('Period (s)')
        ax[2].set_ylabel('Normalised PSD')
        ax[2].legend()
        ax[2].grid()
        plt.tight_layout()
        if i < 10:
            plt.savefig(save_dir +'/simulated_lc_00' + str(i) + '.png')
        else:
            plt.savefig(save_dir +'/simulated_lc_0' + str(i) + '.png')
        plt.clf()







save_dir = '/Users/laura/Documents/QPPs/July_event/july_goddard/detrend_test/smooth_lcs'
def real_data_check():
    smooth_factor = np.arange(11, 301, 6)
    smooth_factor_sg = np.arange(11, 301, 6)
    a, b = fourier_ana(sin_com)
    data_cube = []
    data_cube2 = []
    for i in range(len(smooth_factor)):
        data_cube.append(test_lc - savgol_filter(test_lc, smooth_factor[i]*2+1, 4))
        data_cube2.append(test_lc - smooth(test_lc, smooth_factor[i]))
    for i in range(len(data_cube)):
        f, PSD = fourier_ana(data_cube[i])
        f1, PSD1 = fourier_ana(data_cube2[i])
        fig, ax = plt.subplots(3, figsize = (7,8) )
        ax[0].plot(data_cube[i]/np.max(data_cube[i]), label = 'Savgol window = ' + str(smooth_factor[i]*2+1) + ' s', color = 'g')

        ax[0].grid()
        ax[0].legend()
        ax[0].set_ylabel('Detrended Flux')

        ax[1].plot(data_cube2[i]/np.max(data_cube2[i]), label = 'Smooth window = '+str(smooth_factor[i]), color = 'r')
        ax[1].grid()
        ax[1].legend()
        ax[1].set_ylabel('Detrended Flux')

        ax[2].plot(1./f, PSD/np.max(PSD),label = 'FT', color = 'g')
        ax[2].plot(1./f1, PSD1/np.max(PSD1),label = 'FT', color = 'r')
        
        ax[2].set_xlabel('Period (s)')
        ax[2].set_ylabel('Normalised PSD')
        ax[2].grid()
        plt.tight_layout()
        if i < 10:
            plt.savefig(save_dir +'/both_lc_00' + str(i) + '.png')
        else:
            plt.savefig(save_dir +'/both_lc_0' + str(i) + '.png')
        plt.clf()








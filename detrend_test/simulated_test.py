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
red_noise = np.real(r_n.sample())

#increaing period sin wave
period = np.linspace(10, 50, len(x))
freq = 1./period
amp = np.linspace(1, 0.3, len(x))
amp = np.logspace(3, 0.1, len(x))
amp = gaussian(x, 2000, 800) + 1
sin_com = []
for i in range(len(freq)):
    sin_com.append(amp[i]*np.sin(2*np.pi*freq[i]*x_ar[i]))
sin_com = np.array(sin_com)/1000

ts = SimulatedFlare(flare_background, white_noise, red_noise, sin_com, n, dt).flare()

from wavelet_test import wave_fn
wave_fn(a, 1, cmapp = 'viridis')



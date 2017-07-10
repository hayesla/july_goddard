import numpy as np
import matplotlib.pyplot as plt
from sunpy import lightcurve as lc
import scipy
import matplotlib.dates as dates
from pandas import*
from sunpy.instr import lyra
from astropy.io import fits
import datetime
smooth = scipy.ndimage.filters.uniform_filter 
plt.ion()
from scipy import fftpack
from wave_fn import wave_fn
from scipy.signal import savgol_filter


date = '2012-07-19'
t_start = '2012-07-19 04:15:00'
t_end = '2012-07-19 09:27:00'

#functions#
def normalise(x): #Function to normalise data
	return (x-x.min())/(x.max() - x.min())


def fourier_ana(x):
    N = len(x)
    dt = 2.0
    df = 1./(N*dt)
    PSD = abs(dt*fftpack.fft(x)[:N/2])**2
    f = df*np.arange(N/2)
    f = f[20:]
    PSD = PSD[20:]
    plt.plot(1./f, PSD/np.max(PSD))
    return f, PSD

###-------------DATA SETUP-------------###

#------# 
# GOES #
#------#

goes_lc = lc.GOESLightCurve.create(t_start, t_end)
goeslong = goes_lc.data['xrsb']
goesshort = goes_lc.data['xrsa']

#-----#
# AIA #
#-----#
aia_full = np.loadtxt('aia_131_24s_final.txt', dtype = 'str')
aia_times = []
aia_data = []
for i in range(len(aia_full)):
	att = datetime.datetime.strptime(aia_full[:,0][i] + ' ' + aia_full[:,1][i], '%Y-%m-%d %H:%M:%S.%f')
	aia_times.append(att)
	aia_data.append(float(aia_full[:,2][i]))

aia_lc = Series(aia_data, index = aia_times)

#----------#
# ESP DATA #
#----------#


data_table, header_table = fits.getdata('esp_L1_2012201_005.fit', 1, header = True)

year = data_table.field('year')
doy = data_table.field('doy')
hour = data_table.field('hour')
minute = data_table.field('minute')
sec = data_table.field('sec')

time_frame = []
for i in range(len(year)):
	a = datetime.datetime(2012, 1,1) +datetime.timedelta(days = int(doy[i])-1, hours = int(hour[i]), minutes = int(minute[i]), seconds=float(sec[i]))
	time_frame.append(a)


sxr = data_table.field('QD')
channel18 = data_table.field('CH_18')
channel26 = data_table.field('CH_26')
channel30 = data_table.field('CH_30')
channel36 = data_table.field('CH_36')

sxr = Series(sxr, index = time_frame)
sxr = sxr.truncate(t_start, t_end)
sxr = sxr.resample('2s', how = 'mean')

#####---------- TESTS ---------------####
#				     ####	
#------------------------------------####
test_g = []
x = np.arange(4, 100, 2)
for i in x:
	test_g.append(goeslong.resample(str(i)+'s', how = 'mean'))

new_g = []
for i in range(len(test_g)):
	new_g.append(Series(normalise(np.gradient(test_g[i])), index = test_g[i].index))


def test_timescale(new_g, title = 'Test'):
	detrend = 320.
	dtt = (new_g.index[2] - new_g.index[1]).total_seconds()
	smooth_value = detrend/dtt

	gl_1 = new_g - smooth(new_g, smooth_value)
	gl_2 = smooth(np.abs(gl_1), 2*smooth_value)
	sst = gl_1/gl_2
	time_frame = []
	for i in range(len(new_g)):
		a = datetime.datetime.strptime(str(new_g.index[i]), '%Y-%m-%d %H:%M:%S')
		time_frame.append(a)
	time_frame = np.array(time_frame)
	wave_fn(sst, time_frame, dtt, title)
	
	

### Savitzky-Golay filter ###

box = 201
poly = 3

gl_s = savgol_filter(goeslong, box, poly)

esp_s = savgol_filter(sxr, box, poly)


fig, ax = plt.subplots(2, sharex = True)
ax[0].plot(goeslong, label = 'GOES 1-8 $\mathrm{\AA}$')
ax[0].xaxis.set_major_locator(dates.MinuteLocator(interval =30))
ax[0].xaxis.grid(True, which="major")	
ax[0].xaxis.set_major_formatter(dates.DateFormatter('%H.%M'))
ax[0].legend()
ax[0].set_ylabel('Flux')


ax[1].plot(goeslong - gl_s, label = 'Smoothed by SG filter('+str(box)+','+str(poly)+')')
ax[1].xaxis.set_major_locator(dates.MinuteLocator(interval =30))
ax[1].xaxis.grid(True, which="major")	
ax[1].xaxis.set_major_formatter(dates.DateFormatter('%H.%M'))
ax[1].legend()
ax[1].set_ylabel('Detrended Flux')
ax[1].set_xlabel('Start time ' + str(goeslong.index[0])[0:16] + ' UT')





'''test_final = []
for i in range(1,len(new_g)):
	ts = new_g[i].reindex(new_g[0].index).interpolate()
	test_final.append(ts)

'''

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

date = '2012-07-19'
t_start = '2012-07-19 04:15:00'
t_end = '2012-07-19 09:27:00'

ts = datetime.datetime.strptime(t_start, '%Y-%m-%d %H:%M:%S')
te = datetime.datetime.strptime(t_end, '%Y-%m-%d %H:%M:%S')

#functions#
def normalise(x): #Function to normalise data
	return (x-x.min())/(x.max() - x.min())


#GOES#

goes_lc = lc.GOESLightCurve.create(t_start, t_end)
goeslong = goes_lc.data['xrsb']
goesshort = goes_lc.data['xrsa']

goesl_dif = Series(smooth(np.gradient(goeslong),5), index = goeslong.index)
goess_dif = Series(smooth(np.gradient(goesshort),5), index = goesshort.index)

gl_detrend = goeslong - smooth(goeslong, 41)
gs_detrend = goesshort - smooth(goesshort, 41)

#LYRA#

lyra_full = lc.LYRALightCurve.create(date)
lyra_lc = lyra_full.truncate(t_start, t_end)
lyra_lc_rem = lyra.remove_lytaf_events_from_lightcurve(lyra_lc, ['LAR'], True, '/Users/laura/Documents/QPPs/GI_lyra')[0]

alchaa3 = lyra_lc_rem.data['CHANNEL3']
zrchaa4 = lyra_lc_rem.data['CHANNEL4']

zrcha4 = zrchaa4.resample('2s', how = 'mean')

nanlar = []
for i in range(len(zrcha4)):
	if np.isnan(zrcha4[i]) == True:
        	nanlar.append(i)

zr4 = lyra_lc.data['CHANNEL4']
zr4 = zr4.resample('2s', how = 'mean')
dzr4 = smooth(np.gradient(zr4),7)
dzr4[nanlar] = np.nan
dzr4 = Series(dzr4, index = zr4.index)


#EVE/ESP#

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

#ESP DATA#
sxr = data_table.field('QD')
channel18 = data_table.field('CH_18')
channel26 = data_table.field('CH_26')
channel30 = data_table.field('CH_30')
channel36 = data_table.field('CH_36')

sxr = Series(sxr, index = time_frame)
sxr = sxr.truncate(ts, te)
sxr = sxr.resample('2s', how = 'mean')


#######plotting#########

fig, ax1 = plt.subplots(2, sharex = True)

ax1[0].plot(goeslong.index.to_pydatetime(), normalise(goeslong), label = 'GOES 1-8$\AA$')
ax1[0].plot(goesshort.index.to_pydatetime(), normalise(goesshort), label = 'GOES 0.5-4$\AA$')
#ax1[0].plot(zrcha4.index, normalise(zrcha4), label = 'LYRA Zr')
ax1[0].plot(sxr.index, normalise(sxr), label = 'ESP 0.1-7nm')
#ax1[0].plot(nobe_corr.index, normalise(nobe_corr), label = 'NoRH 17GHz')
ax1[0].set_ylabel('Irradiance (Normalised)', color = 'black', fontsize = 10)

ax1[0].legend(loc = 'upper right', fontsize = 10)
ax1[0].set_ylim(0,1.1)
ax1[0].xaxis.set_major_locator(dates.MinuteLocator(interval = 20))
ax1[0].xaxis.set_major_formatter(dates.DateFormatter('%H.%M'))
ax1[0].xaxis.grid(True, which="major")
ax1[0].set_title('Lightcurves and Derivative Lightcurves', fontsize = 10)


ax1[1].plot(goeslong.index, normalise(goeslong - smooth(goeslong, 100)), label = 'GOES 1-8$\AA$')
#ax1[1].plot(goesshort.index.to_pydatetime(), normalise(goesshort - smooth(goesshort, 100)), label = 'GOES 0.5-4$\AA$')
ax1[1].plot(sxr.index.to_pydatetime(), smooth(normalise(sxr - smooth(sxr, 100)),5)+0.1, label = 'ESP 1-70 $\AA$', color = 'r')




ax1[1].set_ylabel('Detrended Irradiences', color = 'black', fontsize = 10)
ax1[1].xaxis.set_major_locator(dates.MinuteLocator(interval = 20))
ax1[1].xaxis.set_major_formatter(dates.DateFormatter('%H.%M'))
ax1[1].xaxis.grid(True, which="major")
aa = goeslong.index[0]
aa = str(aa)
aa = aa[0:19]
#ax1[1].set_xlabel('Start time:' + aa + 'UT', fontsize = 20)

ax1[1].set_title('')
ax1[1].legend(loc = 'upper right')

'''ax1[2].plot(sxr.index.to_pydatetime(), (sxr - smooth(sxr, 100)), label = 'ESP 0.1-7nm')
ax1[2].plot(dd_d.index, (dd_d), label = 'LYRA Zr', color = 'red')
ax1[2].set_ylabel('Detrended Irradiences', color = 'black', fontsize = 20)
ax1[2].set_ylim(-1e-5, 2e-5)
ax1[2].xaxis.set_major_locator(dates.MinuteLocator(interval = 20))
ax1[2].xaxis.set_major_formatter(dates.DateFormatter('%H.%M'))
ax1[2].xaxis.grid(True, which="major")
aa = goeslong.index[0]
aa = str(aa)
aa = aa[0:19]
ax1[2].set_xlabel('Start time:' + aa + 'UT', fontsize = 20)

ax1[2].set_title('')
ax1[2].legend(loc = 'upper left')'''

'''
import matplotlib as mpl
label_size = 15
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 

fig, ax1 = plt.subplots(1)
ax1.plot(rhessi3_6.index.to_pydatetime(), rhessi3_6 ,label = 'RHESSI 3-6 keV')
ax1.plot(rhessi6_12, label = 'RHESSI 6-12 keV')
ax1.plot(rhessi12_25, label = 'RHESSI 12-25 keV')
ax1.plot(rhessi25_50, label = 'RHESSI 25-50 keV')
ax1.set_yscale('log')
ax1.set_ylim(4, 5000)
ax1.set_ylabel('Count Rate $\mathrm{s^{-1} detector^{-1}}$', color = 'black', fontsize = 20)
ax1.legend(loc = 'upper right', fontsize = 20)
ax1.xaxis.set_major_locator(dates.MinuteLocator(interval =20))
ax1.xaxis.set_major_formatter(dates.DateFormatter('%H.%M'))
ax1.xaxis.grid(True, which="major")
ax1.set_title('RHESSI Lightcurves', fontsize = 20)


aa = str(aa)
aa = aa[0:19]
ax1.set_xlabel('Start time:' + aa + ' UT', fontsize = 20)
plt.tight_layout()


'''


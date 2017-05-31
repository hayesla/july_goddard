import matplotlib.pyplot as plt
import numpy as np
from pandas import Series
from sunpy import map
from sunpy import lightcurve as lc
from scipy import ndimage
smooth = ndimage.filters.uniform_filter
import seaborn as sns
sns.set_style('ticks',{'xtick.direction':'in','ytick.direction':'in'})
sns.set_context('paper', font_scale = 1.5)
import matplotlib.dates as dates

file_dir  = '/Users/laura/Documents/QPPs/July_event/july19/all_aia_fits_24s/*131_.fts'

all_maps = map.Map(file_dir)


def aia_lightcurve():
    aia_data = []
    aia_times = []
    for i in range(len(all_maps)):
        aia_data.append(np.sum(all_maps[i].data))
        aia_times.append(all_maps[i].date)
    aia_lc = Series(aia_data, index = aia_times)
    return aia_lc


def plot_max_value(mapy):
	a = np.where(mapy.data == np.max(mapy.data))
	x_val = a[1][0]
	y_val = a[0][0]
	plt.plot([x_val], [y_val], marker = '.', ms = 10)
	plt.imshow(mapy.data, origin = 'lower')

max_vals = []
for i in range(len(all_maps)):
	a = np.where(all_maps[i].data == np.max(all_maps[i].data))
	x_val = a[1][0]
	y_val = a[0][0]
	max_vals.append([x_val, y_val])


def read_data(filee):
    file_data = np.loadtxt(filee, dtype = 'str')
    dates = []
    x0 = []
    y0 = []
    for i in range(len(file_data)):
        dates.append(datetime.datetime.strptime(file_data[:,0][i] + ' ' + file_data[:,1][i], '%d-%b-%Y %H:%M:%S.%f'))
        x0.append(float(file_data[:,2][i]))
        y0.append(float(file_data[:,3][i]))
    return dates, x0, y0


def plot_together(i):
    
    g = lc.GOESLightCurve.create(all_maps[0].date, all_maps[-1].date)
    gl = g.data['xrsb']
    
	fig = plt.figure(figsize = (12,6))
	ax1 = fig.add_subplot(1,2,1)
	plt.grid(False)
	all_maps[i].plot(vmin = 0, vmax = 16383)
	ax1.set_title('AIA 131 $\mathrm{\AA}$ ' + str(all_maps[i].date)[11:19])

	ax2 = fig.add_subplot(1,2,2)
	ax2.plot(gl.index, normalise(gl))
	ax2.plot(gl.index, normalise(gl - smooth(gl,100)))
	
	ax2.axvline(all_maps[i].date, ls = '--', color = 'k')
	plt.title('Soft X-ray GOES 1-8 $\mathrm{\AA}$')
	ax2.xaxis.set_major_formatter(dates.DateFormatter('%H.%M'))
	ax2.xaxis.set_major_locator(dates.MinuteLocator(interval =30))
	ax2.set_ylabel('Normalised Flux')
	ax2.set_xlabel('Start time 19-July-2012 04:15')


#functions#
def normalise(x): #Function to normalise data
	return (x-x.min())/(x.max() - x.min())


def save_plot_together():
    for i in range(100, len(all_maps)):
        plot_together(i)
        plt.savefig('/Users/laura/Documents/QPPs/July_event/july_goddard/plot_w_lc/aia_'+str(i) + '.png')
        plt.clf()



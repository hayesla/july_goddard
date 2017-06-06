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
from copy import deepcopy
import astropy.units as u

file_dir  = '/Users/laura/Documents/QPPs/July_event/july19/all_aia_fits_24s/*131_.fts'
#rhessi_centroids = '/Users/laura/Documents/QPPs/July_event/july_goddard/96_fwd_fit.csv'
#r_t = read_csv(rhessi_centroids)
#r_times = []
#for i in range(len(r_t)):
#    r_times.append(datetime.datetime.strptime(list(r_t['date'])[i], '%d-%b-%y %H:%M:%S.%f'))

#r_data = list(r_t['data'])


all_maps = map.Map(file_dir)


def aia_lightcurve(map_list):
    aia_data = []
    aia_times = []
    for i in range(len(map_list)):
        aia_data.append(np.sum(map_list[i].data))
        aia_times.append(map_list[i].date)
    aia_lc = Series(aia_data, index = aia_times)
    return aia_lc


def plot_max_value(mapy):
	a = np.where(mapy.data == np.max(mapy.data))
	x_val = a[1][0]
	y_val = a[0][0]
	plt.plot([x_val], [y_val], marker = '.', ms = 10, color = 'k')
	plt.imshow(mapy.data, origin = 'lower', cmap = 'viridis')

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


def read_data_fwdfit(filee):
    file_data = np.loadtxt(filee, dtype = 'str')
    start_time = []
    end_time = []
    flux = []
    x0 = []
    y0 = []
    fwhm0 = []
    for i in range(len(file_data)):
        start_time.append(datetime.datetime.strptime(file_data[:,0][i] + ' ' + file_data[:,1][i], '%Y/%m/%d %H:%M:%S.%f'))
        end_time.append(datetime.datetime.strptime(file_data[:,0][i] + ' ' + file_data[:,2][i], '%Y/%m/%d %H:%M:%S.%f'))
        flux.append(float(file_data[:,3][i]))
        x0.append(float(file_data[:,4][i]))
        y0.append(float(file_data[:,5][i]))
        fwhm0.append(float(file_data[:,6][i]))
    return start_time,end_time,flux, x0, y0, fwhm0



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

def plot_contours():
    test_map = deepcopy(all_maps[0])
    test_map.data = np.zeros(test_map.data.shape)
    test_map.plot(cmap = 'Greys')
    plt.contour(all_maps[320].data, cmap = 'Reds')
    test_map.draw_limb(color =  'k')


#functions#
def normalise(x): #Function to normalise data
	return (x-x.min())/(x.max() - x.min())


def save_plot_together():
    for i in range(100, len(all_maps)):
        plot_together(i)
        plt.savefig('/Users/laura/Documents/QPPs/July_event/july_goddard/plot_w_lc/aia_'+str(i) + '.png')
        plt.clf()



def calc_saturated_lightcurve(new_aia):
    #new_aia = all_maps[320:412]
    sat_lc = []
    new_lc = []
    na_lc = []
    new_aia_index = []
    for i in range(len(new_aia)):
        ma = np.where(new_aia[i].data == np.max(new_aia[i].data))
        na = np.where(new_aia[i].data != np.max(new_aia[i].data))
        val = np.sum(new_aia[i].data[ma[0], ma[1]])
        nal = np.sum(new_aia[i].data[na[0], na[1]])
        sat_lc.append(val)
        na_lc.append(nal)
        new_lc.append(np.sum(new_aia[i].data))
        new_aia_index.append(new_aia[i].date)
    saturated_lc = Series(sat_lc, index = new_aia_index)
    unsaturated_lc = Series(na_lc, index = new_aia_index)
    new_aia_lc = Series(new_lc, index = new_aia_index)
    return saturated_lc, new_aia_lc, unsaturated_lc



def plot_sat():
    new_aia = all_maps[320:412]
    for i in range(0, 10):
        ma = np.where(new_aia[i].data == np.max(new_aia[i].data))
        #fig, ax = plt.subplots(figsize = (6,6))
        plt.grid(False)
        new_aia[i].plot()
        plt.plot(ma[1], ma[0], marker = '.', ls = ' ', color = 'red')
        new_aia[i].plot()
        plt.savefig('/Users/laura/Documents/QPPs/July_event/july_goddard/aia_tests/sat_figs/'+str(i) + '.png')

        plt.clf()


def make_subplots():
    new_aia = all_maps[320:412]
    #new_aia = all_maps
    lengthx = 116 * u.arcsec
    lengthy = 82 * u.arcsec
    x0 = 1016 * u.arcsec
    y0 = -235 * u.arcsec
    new_subs = []
    for i in range(len(new_aia)):
        ssub = submap = new_aia[i].submap(u.Quantity([x0 - lengthx, x0 + lengthx]),
           u.Quantity([y0 - lengthy, y0 + lengthy]))
        new_subs.append(ssub)


def fourier_ana(x):
    N = len(x)
    dt = 24.0
    df = 1./(N*dt)
    PSD = abs(dt*fftpack.fft(x)[:N/2])**2
    f = df*np.arange(N/2)
    #plt.plot(1./f, PSD/np.max(PSD))
    return f, PSD

from scipy import fftpack
def make_freq_map(map_list):
    freq_map = deepcopy(map_list[0])
    full_list = []
    for i in range(len(map_list[0].data[0])):
        for j in range(len(map_list[0].data[:,0])):
            temp_lc = []
            for k in range(len(map_list)):
                temp_lc.append(map_list[k].data[j][i])
        
        
            full_list.append(temp_lc)
        
            temp_test = temp_lc - savgol_filter(temp_lc, 13, 3)
            temp_test = smooth(temp_test,3)
            f, psd = fourier_ana(temp_test)
            p_val = 1./f[np.where(psd == np.max(psd))[0][0]]
            freq_map.data[j][i] = p_val

    return freq_map, full_list

def lala():
    test = np.arange(0, 274*386, 386)
    
    for i in range(len(test)):
        if test[i] != test[len(test)-1]:
            testy.append(full_list[test[i]:test[i+1]])


def chippy_choppy():
    
    #coords found from ginput
    coords = [(48.510204081632622, 78.357606679035257),
    (68.844155844155807, 72.638682745825605),
    (89.813543599257855, 76.451298701298697),
    (115.86641929499069, 82.170222634508349),
    (141.91929499072353, 94.243506493506487),
    (162.25324675324674, 107.58766233766232),
    (173.69109461966602, 127.92161410018551),
    (174.96196660482374, 150.16187384044525),
    (163.52411873840444, 172.40213358070497),
    (142.55473098330239, 190.82977736549162),
    (117.13729128014839, 200.99675324675323),
    (89.178107606679021, 204.17393320964746),
    (68.844155844155807, 207.98654916512055)]

    #make into array for indexing
    coord_array = np.array(coords, dtype = 'int')

    xx = np.linspace(new_subs[0].xrange[0], new_subs[0].xrange[1], 386)[coord_array[:,0]]
    yy = np.linspace(new_subs[0].yrange[0], new_subs[0].yrange[1], 274)[coord_array[:,1]]
    
    
    lc_list = []
    for k in range(len(xx)):
        
        subys = []
        for i in range(len(new_subs)):
            x0 = xx[k]
            y0 = yy[k]
            length = 5*u.arcsec
            subb = new_subs[i].submap(u.Quantity([x0 - length, x0 + length]), u.Quantity([y0 - length, y0 + length]))
            subys.append(subb)

        lc_list.append(aia_lightcurve(subys))
    
    plot_boxes = False
    if plot_boxes:
        new_subs[0].plot()
        for i in range(len(xx)):
            length = 5*u.arcsec
            length2 = length*2
            bottom_left = u.Quantity([xx[i] - length, yy[i] - length])
            new_subs[0].draw_rectangle(bottom_left, length2, length2, linewidth = 1, color = 'k')

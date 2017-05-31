import numpy as np
import matplotlib.pyplot as plt
from sunpy import lightcurve as lc
from scipy import fftpack
from pandas import Series
import scipy
smooth = scipy.ndimage.filters.uniform_filter
from scipy.signal import savgol_filter


def normalise(x):
    return (x-x.min())/(x.max() - x.min())

def fourier_ana(x, dt):
    N = len(x)
    dt = dt
    df = 1./(N*dt)
    PSD = abs(dt*fftpack.fft(x)[:N/2])**2
    f = df*np.arange(N/2)
    plt.plot(1./f, PSD)
    return f, PSD

def fourier_ana_np(x, dt):
    N = len(x)
    f = np.fft.fftfreq(N, dt)
    fp = f > 0
    P = np.absolute(np.fft.fft(x))**2
    freq = f[fp]
    PSD = P[fp]
    plt.plot(1./freq, PSD)
    return freq, PSD

date = '2012-07-19'
t_start = '2012-07-19 04:15:00'
t_end = '2012-07-19 09:27:00'

goes_lc = lc.GOESLightCurve.create(t_start, t_end)
goeslong = goes_lc.data['xrsb']
goesshort = goes_lc.data['xrsa']

class PeriodTest:
    def __init__(self, timeseries):

        '''
        Takes a pandas timeseries and calculates the FT
        '''
            
            
        self.timeseries = timeseries
        
        self.n = len(timeseries)
        
        self.dt = (timeseries.index[1] - timeseries.index[0]).total_seconds()

        self.freq = np.fft.fftfreq(self.n, self.dt)
        
        self.pos_freqs = self.freq[self.freq > 0]
        
        self.PSD = np.absolute(np.fft.fft(timeseries))**2
        
        self.PSD_pos = self.PSD[self.freq>0]


    def plot_fourier(self, logy = False):
        if logy:
            plt.loglog(self.pos_freqs, self.PSD_pos)
            plt.xlabel('Frequency (s)')
            plt.ylabel('PSD')
        else:
            plt.subplot(2,1,1)
            plt.plot(self.timeseries)
            plt.grid()
            plt.subplot(2,1,2)
            plt.plot(1./self.pos_freqs, self.PSD_pos/np.max(self.PSD_pos))
            plt.xlabel('Frequency (s)')
            plt.ylabel('PSD')
            plt.grid()



def difference_test():
    test_g = []
    x = np.arange(4, 100, 2)
    for i in x:
        test_g.append(goeslong.resample(str(i)+'s', how = 'mean'))

    new_g = []
    for i in range(len(test_g)):
        new_g.append(Series(normalise(np.gradient(test_g[i])), index = test_g[i].index))

    nn = []
    for i in range(len(new_g)):
        n_n = len(new_g[i])/23
        if n_n % 2 == 0:
            nn.append(n_n+1)
        else:
            nn.append(n_n)

    ts = []
    for i in range(len(new_g)):
        testy = new_g[i] - savgol_filter(new_g[i], nn[i] ,3)
        ts.append(PeriodTest(testy))

    return ts



def test_chops():
    x = np.arange(0, len(goeslong),500)
    new_lc = []
    test_lc = goeslong - savgol_filter(goeslong, 201, 3)
    for i in range(len(x)-1):
        new_lc.append(test_lc[x[i]:x[i+1]])


    ts = []
    for i in range(len(new_lc)):
        ts.append(PeriodTest(new_lc[i]))

    return ts


def ploty(x, i):
    plt.plot(x[0].timeseries.index, x[i].timeseries/np.max(x[i].timeseries))



import numpy as np
from timeseries_simulation import TimeSeriesFromModelSpectrum
from rnspectralmodels4 import power_law


class SimulatedFlare:
    def __init__(self, background, white_noise, red_noise, periodic_signal, n, dt):
        '''
        create a synthetic solar flare with properties of differnet noise 
        and an underlying periodic component 
        parameters:

        background: polynormial or gaussian; np.array
        white_noise: np.array, normally distributed noise
        red_noise: power_law distribution 
        periodic_signal: the underlying periodic component
        n: number of points
        dt: sampling rate. 
        '''

        self.background = background

        self.white_noise = white_noise

        self.red_noise = red_noise

        self.periodic_signal = periodic_signal

        self.n = n

        self.dt = dt

    def flare(self):
        return self.background + self.white_noise + self.red_noise + self.periodic_signal


    def plot_flare(self):
        fig, ax = plt.subplots()
        ax.plot(self.flare(), label = 'Simulated Flare')
        ax.set_xlabel('X values')
        ax.set_ylabel('Normalised Simulated Flux')
        ax.set_title('Background + white noise + red noise + periodic signal')
        ax.set_ylim(np.min(self.flare()) - 0.05, np.max(self.flare())+0.05)
        plt.tight_layout()

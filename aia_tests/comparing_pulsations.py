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
from scipy.signal import savgol_filter
from aia_tests import aia_lightcurve, calc_saturated_lightcurve

all_maps = map.Map(file_dir)

g = lc.GOESLightCurve.create(all_maps[0].date, all_maps[-1].date)
gl = g.data['xrsb']

aia_lc = aia_lightcurve()

gl_test = gl - savgol_filter(gl, 201, 3)
aia_test = aia_lc - savgol_filter(aia_lc, 13, 3)

sat, new, unsat = calc_saturated_lightcurve()
sat_test = sat - savgol_filter(sat, 13, 3)
new_test = new - savgol_filter(new, 13, 3)
unsat_test = unsat - savgol_filter(unsat, 13, 3)


fig, ax = plt.subplots(2, sharex = True)
ax[0].plot(aia_lc/np.max(aia_lc), label = 'AIA 131 A')
ax[0].plot(gl/np.max(gl), color = 'r', label = 'GOES 1-8 A')
ax[0].legend()


ax[1].plot(aia_test/np.max(aia_test), label = 'AIA 131 A')
ax[1].plot(gl_test/np.max(gl_test), label = 'GOES 1-8 A', color = 'r')

ax[1].plot(new_test/np.max(new_test))

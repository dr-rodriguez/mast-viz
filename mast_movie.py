# Read the HDF5 stores and generate frames for a movie
import matplotlib
matplotlib.use('Qt5Agg')  # avoids crashing MacOS Mojave
import numpy as np
import pandas as pd
import healpy as hp
import copy
from matplotlib import cm
from astropy.time import Time
import matplotlib.pyplot as plt
plt.interactive(False)

# GALEX
rootname = 'galex'
df = pd.read_hdf('data/galex.h5', 'data')
ptab = pd.read_hdf('data/galex.h5', 'ptab')
base_map = hp.read_map('data/galex_map.fits')
exp_map = np.zeros(len(base_map))

# Min and Maximum time range
tstep = 7.  # days
t0 = np.min(df['t_min'])
tf = np.max(df['t_min'])

# Weekly bins
week_bin = np.trunc((df['t_min'] - t0) / tstep)
df['week_bin'] = week_bin
weeks = df.groupby('week_bin')

MIN, MAX = [], []

# Plotting setup
plt.style.use('dark_background')
cmap = copy.copy(cm.get_cmap('cividis'))
cmap.set_bad('xkcd:charcoal')
cmap.set_under('navy')
plt.rcParams.update({'font.size': 15})
lon = np.arange(360)
lat = np.zeros(360)

for w in weeks.groups.keys():
    week_data = weeks.get_group(w)
    print(w, len(week_data))
    title = ''
    for _, row in week_data.iterrows():
        # Get the time for the plot title
        if title == '':
            t1 = row['t_min']
            tobj = Time(t1, format='mjd')
            title = tobj.datetime.strftime('%Y')

        # Get the pixels for that observation and add them to the map
        pix = ptab[ptab['obs_id'] == row['obs_id']]
        pix = pix.iloc[0]  # in case there are duplicate rows
        exp_map[pix['ind']] = exp_map[pix['ind']] + row['t_exptime']

    smap = exp_map.copy()
    smap = np.ma.masked_where(smap == 0, smap)
    smap.set_fill_value(np.nan)

    MIN_ = smap.min()
    MAX_ = smap.max()
    MIN.append(MIN_)
    MAX.append(MAX_)

    hp.mollview(np.log10(smap), cmap=cmap, rot=-80,
                flip='geo', coord='C', cbar=False, notext=True,
                norm='linear', xsize=1000, title=title)
    hp.projplot(lon, lat, 'r', lonlat=True, coord='G')
    hp.graticule(dpar=45., dmer=30., coord='C')
    pngfile1 = 'movie/frame/' + rootname + '_frame' + str(i) + '.png'
    plt.savefig(pngfile1, dpi=300)
    plt.close()

print("min(MIN), max(MAX) = ", min(MIN), max(MAX))

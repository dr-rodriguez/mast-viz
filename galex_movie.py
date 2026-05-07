# Read the HDF5 stores and generate frames for a movie
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd
import healpy as hp
import copy
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap
from astropy.time import Time
import matplotlib.pyplot as plt
plt.interactive(False)

# GALEX
movie_dir = 'movie/galex/'
rootname = 'galex'
DPI=100  # 300 is default
HIGHLIGHTS = True  # flag for flashing new observations as they are added
GALACTIC_LINE = False  # flag for enabling red galactic line
SKYCOLOR = '#003B4D'  # MAST darkest turquoise
# SKYCOLOR = 'midnightblue'

df = pd.read_hdf('data/galex.h5', 'data')
ptab = pd.read_hdf('data/galex.h5', 'ptab')
base_map = hp.read_map('data/galex_map.fits')

# Min and Maximum time range
tstep = 7.  # days
t0 = np.min(df['t_min'])
tf = np.max(df['t_min'])
date0 = df[df['t_min'] == t0]['t_min'].iloc[0]

# Weekly bins
week_bin = np.trunc((df['t_min'] - t0) / tstep)
df['week_bin'] = week_bin
weeks = df.groupby('week_bin')

# Plotting setup
plt.style.use('dark_background')
cmap = copy.copy(cm.get_cmap('cividis'))
# cmap.set_bad('xkcd:charcoal')
cmap.set_bad(SKYCOLOR)
cmap.set_under('k')
# Color "map" for highlights:
# color_array = np.array([[240./256., 199./256., 36./256., 0],  # RGBA, A=0 is transparent
#                         [240./256., 199./256., 36./256., 1]])
# MAST Orange C75109 (199 81 9)
color_array = np.array([[199./256., 81./256., 9./256., 0],  # RGBA, A=0 is transparent
                        [199./256., 81./256., 9./256., 1]])
highlight_cmap = LinearSegmentedColormap.from_list(name='highlight', colors=color_array)
plt.rcParams.update({'font.size': 15})
lon = np.arange(360)
lat = np.zeros(360)

# Some initial setup
exp_map = np.zeros(len(base_map))
MIN, MAX = [], []
time_range = range(0, int(max(week_bin))+1)

# Store some statistics with time
time_stats = []

# Main loop
for i in time_range:
    w = time_range[i]
    area = 0
    obs_counts = 0
    exp_counts = 0
    highlights = np.zeros(len(base_map))  # reset highlights map
    try:
        week_data = weeks.get_group(w)
        print(i, w, len(week_data))
        title = ''
        for _, row in week_data.iterrows():
            # Get the time for the plot title
            if title == '':
                t1 = row['t_min']
                tobj = Time(t1, format='mjd')
                title = tobj.datetime.strftime('%Y-%m')

            # Get the pixels for that observation and add them to the map
            pix = ptab[ptab['obs_id'] == row['obs_id']]

            # Skip missing data (eg, bad footprints)
            if len(pix) == 0:
                continue

            pix = pix.iloc[0]  # in case there are duplicate rows
            exp_map[pix['ind']] = exp_map[pix['ind']] + row['t_exptime']
            highlights[pix['ind']] = exp_map[pix['ind']] + 1  # for highlighting recent additions

            obs_counts += 1
            area += len(pix['ind'])
            exp_counts += row['t_exptime']
    except KeyError:
        # no data for this week
        tobj = Time(date0 + tstep*i, format='mjd')
        title = tobj.datetime.strftime('%Y-%m')
        pass

    # Populate time-based stats
    datum = {'week': w, 'date': tobj, 'area': area, 'obs_counts': obs_counts, 'exp_counts': exp_counts}
    time_stats.append(datum)

    smap = exp_map.copy()
    smap = np.ma.masked_where(smap == 0, smap)
    smap.set_fill_value(np.nan)

    MIN_ = smap.min()
    MAX_ = smap.max()
    MIN.append(MIN_)
    MAX.append(MAX_)

    # Viewing angle, from -180 to 180
    rot_angle = (i % 360) - 180

    # Make the actual plot
    hp.mollview(np.ma.log10(smap), cmap=cmap, rot=rot_angle,
                min=0.0, max=7.,  # exptime limits
                flip='geo', coord='C',
                cbar=False, notext=True,
                bgcolor='black', badcolor=SKYCOLOR,
                norm='linear', xsize=1000, title=title)

    if HIGHLIGHTS:
        # Highlight new observations this week
        hmap = highlights.copy()
        hmap = np.ma.masked_where(hmap == 0, hmap)
        hmap.set_fill_value(np.nan)

        hp.mollview(np.ma.log10(hmap), cmap=highlight_cmap, rot=rot_angle,
                    min=0.0, max=1.0,
                    flip='geo', coord='C',
                    cbar=False, notext=True,
                    bgcolor=[0., 0., 0., 0.], badcolor=[0., 0., 0., 0.],  # fully transparent colors
                    norm='linear', xsize=1000, title=title,
                    reuse_axes=True, fig=1)

    # Galatic line
    if GALACTIC_LINE:
        hp.projplot(lon, lat, 'r', lonlat=True, coord='G')

    # Grid
    hp.graticule(dpar=45., dmer=30., coord='C', color='lightgray')

    # Output file
    pngfile1 = movie_dir + rootname + f'_frame{i:06d}.png'
    plt.savefig(pngfile1, dpi=DPI)
    plt.close()

print("min(MIN), max(MAX) = ", min(MIN), max(MAX))

# Write out time stats
time_df = pd.DataFrame(time_stats)
time_df.to_csv('data/' + rootname + '_time.csv', index=False)

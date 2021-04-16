# Read the HDF5 stores and generate frames for a movie
import matplotlib

matplotlib.use('Qt5Agg')  # avoids crashing MacOS Mojave
import os
import numpy as np
import pandas as pd
import healpy as hp
import copy
from matplotlib import cm
from astropy.time import Time
import matplotlib.pyplot as plt

plt.interactive(False)

# PS1
print('Loading data...')
rootname = 'ps1'
movie_dir = 'movie/ps1/'
df = pd.read_hdf('data/ps1.h5', 'data')
ptab = pd.read_hdf('data/ps1.h5', 'ptab')
base_map = hp.read_map('data/ps1_map.fits')

# Min and Maximum time range
tstep = 7.  # days
t0 = np.min(df['t_min'])
tf = np.max(df['t_min'])
date0 = df[df['t_min'] == t0]['t_min'].iloc[0]

# Weekly bins
print('Generating week bins...')
week_bin = np.trunc((df['t_min'] - t0) / tstep)
df['week_bin'] = week_bin
weeks = df.groupby('week_bin')

# Plotting setup
plt.style.use('dark_background')
cmap = copy.copy(cm.get_cmap('cividis'))
cmap.set_bad('xkcd:charcoal')
cmap.set_under('k')
plt.rcParams.update({'font.size': 15})
lon = np.arange(360)
lat = np.zeros(360)

# Some initial setup
exp_map = np.zeros(len(base_map))
MIN, MAX = [], []
time_range = range(0, int(max(week_bin)) + 1)

# Store some statistics with time
time_stats = []

print('Starting image loop')
print(f'{len(time_range)} steps to process...')
resume = False
w_resume = 0

# Check if this is a continuing run and re-start from there
if os.path.isfile(movie_dir + 'temp_time.csv') and os.path.isfile(movie_dir + 'temp_expmap.npy'):
    print('Resuming from last run')

    # Read temp_time file, convert to a list of dictionaries, and get the maximum week processed
    temp_time = pd.read_csv(movie_dir + 'temp_time.csv')
    time_stats = temp_time.to_dict('records')
    w_resume = temp_time['week'].max()
    
    # Read temp_smap file
    exp_map = np.load(movie_dir + 'temp_expmap.npy')

    resume = True

# Main loop
for i in time_range:
    w = time_range[i]
    
    # Resume after the last processed week
    if resume and w <= w_resume:
        continue
        
    area = 0
    obs_counts = 0
    exp_counts = 0
    try:
        week_data = weeks.get_group(w)
        print(i, w, len(week_data), Time.now())
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

            obs_counts += 1
            area += len(pix['ind'])
            exp_counts += row['t_exptime']
    except KeyError:
        # no data for this week
        tobj = Time(date0 + tstep * i, format='mjd')
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
    hp.mollview(np.log10(smap), cmap=cmap, rot=rot_angle,
                min=0.0, max=7.,  # exptime limits
                flip='geo', coord='C',
                cbar=False, notext=True,
                bgcolor='black', badcolor='midnightblue',
                norm='linear', xsize=1000, title=title)
    hp.projplot(lon, lat, 'r', lonlat=True, coord='G')
    hp.graticule(dpar=45., dmer=30., coord='C', color='lightgray')
    pngfile1 = movie_dir + rootname + f'_frame{i:06d}.png'
    plt.savefig(pngfile1, dpi=300)
    plt.close()

    # Temporary save of time_stats and exp_map to be able to resume processing
    time_df = pd.DataFrame(time_stats)
    time_df.to_csv(movie_dir + 'temp_time.csv', index=False)
    np.save(movie_dir + 'temp_expmap.npy', exp_map)

print("min(MIN), max(MAX) = ", min(MIN), max(MAX))

# Write out time stats
time_df = pd.DataFrame(time_stats)
time_df.to_csv('data/' + rootname + '_time.csv', index=False)

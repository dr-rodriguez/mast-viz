import os
import copy
import warnings
import numpy as np
import pandas as pd
import healpy as hp
from astropy.time import Time
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ruff: noqa: PLC0415
# Suppress healpy warning about existing figures
warnings.filterwarnings("ignore", category=UserWarning, message=".*Ignoring specified arguments.*")

def setup_colormaps(skycolor):
    try:
        from matplotlib import colormaps
        cmap = copy.copy(colormaps['cividis'])
    except Exception:
        from matplotlib import cm
        cmap = copy.copy(cm.get_cmap('cividis'))
    
    cmap.set_bad(skycolor)
    cmap.set_under('k')

    # Color "map" for highlights:
    # MAST Orange C75109 (199 81 9)
    color_array = np.array([[199./256., 81./256., 9./256., 0],  # RGBA, A=0 is transparent
                            [199./256., 81./256., 9./256., 1]])
    highlight_cmap = LinearSegmentedColormap.from_list(name='highlight', colors=color_array)
    return cmap, highlight_cmap


def get_time_bins(df, tstep):
    t0 = np.min(df['t_min'])
    date0 = df[df['t_min'] == t0]['t_min'].iloc[0]
    
    week_bin = np.trunc((df['t_min'] - t0) / tstep)
    df['week_bin'] = week_bin
    weeks = df.groupby('week_bin')
    return weeks, date0, week_bin


def load_resume_state(movie_dir):
    time_stats = []
    w_resume = -1
    exp_map = None
    if os.path.isfile(movie_dir + 'temp_time.csv') and os.path.isfile(movie_dir + 'temp_expmap.npy'):
        print('Resuming from last run')
        temp_time = pd.read_csv(movie_dir + 'temp_time.csv')
        time_stats = temp_time.to_dict('records')
        w_resume = temp_time['week'].max()
        exp_map = np.load(movie_dir + 'temp_expmap.npy')
    return time_stats, w_resume, exp_map


def save_resume_state(movie_dir, time_stats, exp_map):
    time_df = pd.DataFrame(time_stats)
    time_df.to_csv(movie_dir + 'temp_time.csv', index=False)
    np.save(movie_dir + 'temp_expmap.npy', exp_map)


def plot_frame(smap, rot_angle, title, cmap, skycolor, highlights_map, highlight_cmap, 
               highlights, galactic_line, galactic_grid, movie_dir, rootname, i, dpi):
    _ = plt.figure(1, figsize=(8.5, 5.4))
    
    # Use np.ma.log10 for masked arrays or np.log10 if not masked
    if hasattr(smap, 'mask'):
        log_smap = np.ma.log10(smap)
    else:
        log_smap = np.log10(smap)

    hp.mollview(log_smap, cmap=cmap, rot=rot_angle,
                min=0.0, max=7.,  # exptime limits
                flip='geo', coord='C',
                cbar=False, notext=True,
                bgcolor='black', badcolor=skycolor,
                norm='linear', xsize=1000, title=title,
                fig=1)

    if highlights:
        # Highlight new observations this week
        hmap = highlights_map.copy()
        hmap = np.ma.masked_where(hmap == 0, hmap)
        hmap.set_fill_value(np.nan)
        if hasattr(hmap, 'mask'):
            log_hmap = np.ma.log10(hmap)
        else:
            log_hmap = np.log10(hmap)

        hp.mollview(log_hmap, cmap=highlight_cmap, rot=rot_angle,
                    min=0.0, max=1.0,
                    flip='geo', coord='C',
                    cbar=False, notext=True,
                    bgcolor=[0., 0., 0., 0.], badcolor=[0., 0., 0., 0.],  # type: ignore
                    norm='linear', xsize=1000, title=title,
                    reuse_axes=True, fig=1)

    # Galactic line
    if galactic_line:
        lon = np.arange(360)
        lat = np.zeros(360)
        hp.projplot(lon, lat, 'r', lonlat=True, coord='G')

    # Grid
    if galactic_grid:
        hp.graticule(dpar=45., dmer=30., coord='C', color='lightgray')

    # Output file
    pngfile = movie_dir + rootname + f'_frame{i:06d}.png'
    plt.savefig(pngfile, dpi=dpi)
    plt.close()


def generate_movie(df, ptab, base_map, rootname, movie_dir, 
                   tstep=7.0, 
                   highlights=False, 
                   galactic_line=False, 
                   galactic_grid=True, 
                   skycolor='#003B4D', 
                   dpi=100, 
                   resume=True):
    """
    Main entry point for generating movie frames.
    
    Args:
        df (pd.DataFrame): Main dataset with 't_min' and 't_exptime'.
        ptab (pd.DataFrame or dict): Footprint pixel tables. If dict, expected keys are 'obs_collection'.
        base_map (np.ndarray): Healpy base map for size reference.
        rootname (str): Prefix for generated files.
        movie_dir (str): Directory to save frames.
        tstep (float): Time step in days.
        highlights (bool): Whether to flash new observations.
        galactic_line (bool): Whether to show the galactic plane.
        galactic_grid (bool): Whether to show the coordinate grid.
        skycolor (str): Background color for empty pixels.
        dpi (int): DPI for output images.
        resume (bool): Whether to attempt to resume from temporary state.
    """
    # Ensure movie directory exists
    if movie_dir:
        os.makedirs(movie_dir, exist_ok=True)

    plt.style.use('dark_background')
    plt.rcParams.update({'font.size': 15})
    
    cmap, highlight_cmap = setup_colormaps(skycolor)
    weeks, date0, week_bin = get_time_bins(df, tstep)
    
    time_range = range(0, int(max(week_bin)) + 1)
    
    time_stats = []
    w_resume = -1
    exp_map = np.zeros(len(base_map))
    
    if resume:
        t_stats, w_res, e_map = load_resume_state(movie_dir)
        if w_res >= 0:
            time_stats = t_stats
            w_resume = w_res
            if e_map is not None:
                exp_map = e_map

    MIN, MAX = [], []
    tobj = None
    
    print('Starting image loop')
    print(f'{len(time_range)} steps to process...')
    
    for i in time_range:
        w = time_range[i]
        
        # Resume after the last processed week
        if resume and w <= w_resume:
            continue
            
        area = 0
        obs_counts = 0
        exp_counts = 0
        highlights_map = np.zeros(len(base_map))  # reset highlights map
        
        try:
            week_data = weeks.get_group(w)
            print(i, w, len(week_data), Time.now())
            title = ''
            for _, row in week_data.iterrows():
                # Get the time for the plot title
                if title == '':
                    t1 = row['t_min']
                    tobj = Time(t1, format='mjd')
                    title = tobj.strftime('%Y-%m')

                # Get the pixels for that observation
                if isinstance(ptab, dict):
                    coll = row.get('obs_collection')
                    if coll and coll in ptab:
                        pix = ptab[coll][ptab[coll]['obs_id'] == row['obs_id']]
                    else:
                        continue
                else:
                    pix = ptab[ptab['obs_id'] == row['obs_id']]
                    
                # Skip missing data (eg, bad footprints)
                if len(pix) == 0:
                    continue

                pix = pix.iloc[0]  # in case there are duplicate rows
                exp_map[pix['ind']] = exp_map[pix['ind']] + row['t_exptime']
                highlights_map[pix['ind']] = exp_map[pix['ind']] + 1  # for highlighting recent additions

                obs_counts += 1
                area += len(pix['ind'])
                exp_counts += row['t_exptime']
        except KeyError:
            # no data for this week
            tobj = Time(date0 + tstep * i, format='mjd')
            title = tobj.strftime('%Y-%m')

        # Populate time-based stats
        datum = {'week': w, 'date': tobj, 'area': area, 'obs_counts': obs_counts, 'exp_counts': exp_counts}
        time_stats.append(datum)

        smap = exp_map.copy()
        smap = np.ma.masked_where(smap == 0, smap)
        smap.set_fill_value(np.nan)

        MIN.append(smap.min())
        MAX.append(smap.max())

        # Viewing angle, from -180 to 180
        rot_angle = (i % 360) - 180

        # Render the frame
        plot_frame(smap, rot_angle, title, cmap, skycolor, highlights_map, highlight_cmap, 
                   highlights, galactic_line, galactic_grid, movie_dir, rootname, i, dpi)

        if resume:
            save_resume_state(movie_dir, time_stats, exp_map)

    if len(MIN) > 0 and len(MAX) > 0:
        print("min(MIN), max(MAX) = ", min(MIN), max(MAX))

    # Write out time stats
    time_df = pd.DataFrame(time_stats)
    time_df.to_csv('data/' + rootname + '_time.csv', index=False)

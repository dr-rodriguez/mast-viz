import matplotlib
matplotlib.use('Qt5Agg')  # avoids crashing MacOS Mojave
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.table import Table
import time
from matplotlib import cm
import copy
from utils import parse_s_region, get_polygon_cdshealpix
from db_connect import get_db_data
plt.interactive(False)


def save_db_data(df, datafile='data.csv', fmt='fits'):
    if fmt.lower() == 'csv':
        # Output as CSV file
        df.to_csv(datafile, index=False)
    else:
        # Output as FITS file
        t2 = Table.from_pandas(df, index=False)
        t2.write(datafile, overwrite=True)


def read_file_data(datafile='data.fits', fmt='table'):
    # Read the data, either as a csv into pandas or fits into astropy Table
    if datafile.endswith('fits'):
        t = Table.read(datafile)
    else:
        t = Table.read(datafile, format='csv')
        # df = pd.read_csv(datafile)

    # Decode bytes columns
    for col in t.colnames:
        if isinstance(t[col][0], bytes):
            t[col] = [x.decode() for x in t[col]]

    if fmt.lower() == 'table':
        return t
    else:
        return t.to_pandas()


def make_map(df, nside=256, exp_col='t_exptime', verbose=False):
    """
    Make Healpix map

    Parameters
    ----------
    df
        Source dataframe
    nside
        HEALPix resolution
    exp_col
        Exposure column name

    Returns
    -------
    hp_map
        Exposure map by HEALPix
    ptab
        HEALPix pixel values for the provided dataframe
    """

    # number of pixels for that resolution
    npix = hp.nside2npix(nside)
    hp_map = np.zeros(npix)
    resolution = hp.nside2resol(nside, arcmin=True)
    print(f'NSIDE={nside} NPIX={npix} Resolution(arcmin)={resolution}')

    # I believe this is an estimate for how many HP-resolution elements there are in an area,
    # which for us is very variable. I'm using the TESS FFI size of just over 12 degrees
    # However, this will lead to gigantic fits as the way the output is produced is to store
    # zeros for all possible ranges
    # We may want to construct a sparser data structure which we can then read and interpret
    size_deg = 13  # largest area to consider, in degrees
    nind = int(size_deg ** 2 * np.pi / (resolution / 60) ** 2 * size_deg)

    ptab = []
    # ptab = np.zeros(len(df), dtype=[('np', 'i8'), ('ind', '%di8' % nind)])  # table to store healpix ids

    for i, row in df.iterrows():
        if row['s_region'] is not None and row[exp_col] is not None:

            # Parse the footprint
            try:
                coords = parse_s_region(row['s_region'])  # converts CIRCLE to 16-point POLYGON
                ra_list, dec_list = coords['ra'], coords['dec']
            except:
                print('Unable to parse s_region for {}'.format(row['obs_id']))
                continue

            # Generate HEALPix indices for the footprint
            try:
                # CDSHEALPIX package
                # Indices of all pixels inside or intersecting the polygon
                ipix = get_polygon_cdshealpix(ra_list, dec_list, depth=int(np.log2(nside)))
                if verbose: print('CDS', row['obs_id'], ipix)
            except:
                # HEALPY package
                try:
                    # Indices of all pixels inside or intersecting the convex polygon vec (if not convex it will fail)
                    vec = hp.ang2vec(ra_list, dec_list, lonlat=True)
                    ipix = hp.query_polygon(nside, vec, inclusive=True)
                    if verbose: print('HP ', row['obs_id'], ipix)
                except:
                    print('Unable to get HP indeces for {}: {}'.format(row['obs_id'], row['s_region']))
                    continue

            # Populate output lists
            try:
                hp_map[ipix] = hp_map[ipix] + np.float(row[exp_col])  # adding up exposures for that healpix pixel

                p_dict = {'i': i, 'obs_id': row['obs_id'],  # identifiers
                          'np': len(ipix),  # number of pixels for each observation
                          'ind': ipix}  # pixel values for each observation
                ptab.append(p_dict)
                # ptab['np'][i] = len(ipix)  # number of pixels for each observation
                # ptab['ind'][i][0:len(ipix)] = ipix  # pixel values for each observation
            except:
                print('Unable to store indices for {}: {}'.format(row['obs_id'], ipix))
                continue

            # if row['s_region'].startswith('CIRCLE'):
            #     # Use circle approach (center + radius)
            #
            #     vec = hp.ang2vec(ra0, dec0, lonlat=True)  # healpix vector for position
            #     ipix = hp.query_disc(nside, vec, radius)  # indices of all healpix pixels in a radius around vec
            # else:
            #     # Use polygon approach

    pdf = pd.DataFrame(ptab)

    return hp_map, pdf


def output_map(hp_map, outfile='mast_map.fits'):
    # Output the map
    hp.write_map(outfile, hp_map, coord='C', overwrite=True)
    return


def read_map(mapfile):
    return hp.read_map(mapfile)


def make_plot(hp_map, outfile='mast_map.png', title=''):
    # Generate the map

    # pngfile1 = os.path.splitext(outfile)[0] + '_300.png'
    pngfile2 = os.path.splitext(outfile)[0] + '_600.png'

    # Plot options
    plt.style.use('dark_background')
    cmap = copy.copy(cm.get_cmap('cividis'))
    # cmap.set_under('w')
    cmap.set_bad('xkcd:charcoal')
    cmap.set_under('navy')

    plt.rcParams.update({'font.size': 15})
    lon = np.arange(360)
    lat = np.zeros(360)

    # datestr = time.strftime('%B %d, %Y')
    hp.mollview(np.log10(hp_map + 0.1), cmap=cmap,
                min=0.0, max=6.,  # exptime limits
                rot=-80,
                flip='geo', coord='C',
                cbar=False,  # color bar
                notext=True,
                title=title,
                bgcolor='black', badcolor='gray',
                norm='linear',
                xsize=1000)
    hp.projplot(lon, lat, 'r', lonlat=True, coord='G')
    hp.graticule(dpar=45., dmer=30., coord='C', color='white')
    # plt.savefig(pngfile1, dpi=300)
    plt.savefig(pngfile2, dpi=600)

    plt.close()

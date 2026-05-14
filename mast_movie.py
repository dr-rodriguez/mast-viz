# Read the HDF5 stores and generate frames for a movie
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
import movie_utils

plt.interactive(False)

# MAST
print('Loading data...')
rootname = 'mast_v4'
movie_dir = 'movie/mast_v4/'
HIGHLIGHTS = True
GALACTIC_LINE = False
SKYCOLOR = '#003B4D'
tstep = 7.

df1 = pd.read_hdf('data/galex.h5', 'data')
df2 = pd.read_hdf('data/hst.h5', 'data')
df3 = pd.read_hdf('data/tess.h5', 'data')

df = pd.concat([df1, df2, df3])
df = df.reset_index(drop=True)

# Use a dictionary for ptab to optimize collection-based lookup in mast
ptab = {
    'GALEX': pd.read_hdf('data/galex.h5', 'ptab'),
    'HST': pd.read_hdf('data/hst.h5', 'ptab'),
    'TESS': pd.read_hdf('data/tess.h5', 'ptab')
}

base_map = hp.read_map('data/galex_map.fits')

movie_utils.generate_movie(
    df=df,
    ptab=ptab,
    base_map=base_map,
    rootname=rootname,
    movie_dir=movie_dir,
    tstep=tstep,
    highlights=HIGHLIGHTS,
    galactic_line=GALACTIC_LINE,
    galactic_grid=True,
    skycolor=SKYCOLOR,
    dpi=300,
    resume=True
)

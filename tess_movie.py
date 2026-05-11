# Read the HDF5 stores and generate frames for a movie
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import healpy as hp
import movie_utils

# TESS
movie_dir = 'movie/tess/'
rootname = 'tess'
DPI = 100
HIGHLIGHTS = False
GALACTIC_LINE = False
GALACTIC_GRID = True
SKYCOLOR = '#003B4D'
tstep = 7.

print('Loading data...')
df = pd.read_hdf('data/tess.h5', 'data')
ptab = pd.read_hdf('data/tess.h5', 'ptab')
base_map = hp.read_map('data/tess_map.fits')

movie_utils.generate_movie(
    df=df,
    ptab=ptab,
    base_map=base_map,
    rootname=rootname,
    movie_dir=movie_dir,
    tstep=tstep,
    highlights=HIGHLIGHTS,
    galactic_line=GALACTIC_LINE,
    galactic_grid=GALACTIC_GRID,
    skycolor=SKYCOLOR,
    dpi=DPI,
    resume=True
)

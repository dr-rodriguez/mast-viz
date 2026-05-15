# Read the HDF5 stores and generate frames for a movie
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import healpy as hp
from mast_viz.utils import movie_utils

if __name__ == "__main__":
    # HST
    movie_dir = 'movie/hst/'
    rootname = 'hst'
    DPI = 100
    HIGHLIGHTS = True
    GALACTIC_LINE = False
    GALACTIC_GRID = True
    SKYCOLOR = '#003B4D'
    tstep = 7.

    print('Loading data...')
    df = pd.read_hdf('data/hst.h5', 'data')
    ptab = pd.read_hdf('data/hst.h5', 'ptab')
    base_map = hp.read_map('data/hst_map.fits')

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

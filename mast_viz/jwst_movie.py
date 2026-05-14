# Read the HDF5 stores and generate frames for a movie
import matplotlib
matplotlib.use('Qt5Agg')  # avoids crashing MacOS Mojave
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
from .utils import movie_utils

if __name__ == "__main__":
    plt.interactive(False)

    # JWST
    print('Loading data...')
    rootname = 'jwst'
    movie_dir = 'movie/jwst/'
    tstep = 1.

    df = pd.read_hdf('data/jwst.h5', 'data')
    ptab = pd.read_hdf('data/jwst.h5', 'ptab')
    base_map = hp.read_map('data/jwst_map.fits')

    movie_utils.generate_movie(
        df=df,
        ptab=ptab,
        base_map=base_map,
        rootname=rootname,
        movie_dir=movie_dir,
        tstep=tstep,
        highlights=False,
        galactic_line=True,
        galactic_grid=True,
        skycolor='midnightblue',
        dpi=300,
        resume=True
    )

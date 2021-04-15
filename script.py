# Script to test things out
from mast_plot import *

# MOC order = HEALPix resolution order, the logarithm base 2 of nside
# MOC order 8 -> NSIDES 256
# NOTE: MOC uses NESTED and HP by default uses RING
# nest2ring(nside, ipix)
# Convert pixel number from NESTED ordering to RING ordering.
# A quick investigation suggests that
# moc.serialize()[1].data is a table with column UNIQ which is the unique
# healpix integer for the max_order specified
# Correction: it maybe that it's a combination of ORDER and HP as provided by:
# https://github.com/cds-astro/mocpy/blob/master/mocpy/utils.py
# which means I need to untangle that?
# I found some very useful information with the core CDSHEALPIX package.
# That may be the better tool to use rather than MOCs


# region HST
# For HST I just get imaging data and exclude HAP
# Get from DB
df = get_db_data(mission='HST', limit=None,
                  constraints="AND dataproduct_type='image' AND provenance_name!='HAP' ",
                  server='mastdbtest', db='CAOMv230ProdSync')

# Using HDF5 to handle the complex data storage without using CSV/FITS (except for the map)
store = pd.HDFStore("data/hst.h5")
store['data'] = df

# Get from H5 rather than DB
df = pd.read_hdf("data/hst.h5", "data")

# Generate the map and save it
hp_map, ptab = make_map(df)
store['ptab'] = ptab
output_map(hp_map, outfile='data/hst_map.fits')

# Read from the map file
hp_map = read_map(mapfile='data/hst_map.fits')
store.close()

# Make the plot and save it
make_plot(hp_map, outfile='image/mast_hst_map.png', title='HST')
# endregion

# region TESS
df = get_db_data(mission='TESS', limit=None,
                 constraints="AND dataproduct_type='image' AND sequence_number < 32 ",
                 server='mastdbtest', db='CAOMv230ProdSync')

# Using HDF5 to handle the complex data storage without using CSV/FITS (except for the map)
store = pd.HDFStore("data/tess.h5")
store['data'] = df

# Get from the HDF5 rather than DB
df = pd.read_hdf("data/tess.h5", "data")

hp_map, ptab = make_map(df)

store['ptab'] = ptab
output_map(hp_map, outfile='data/tess_map.fits')  # in FITS format

hp_map = read_map(mapfile='data/tess_map.fits')
store.close()

make_plot(hp_map, outfile='image/mast_tess_map.png', title='TESS')
# endregion

# region GALEX
df = get_db_data(mission='GALEX', limit=None,
                 constraints="AND dataproduct_type='image' ",
                 server='mastdbtest', db='CAOMv230ProdSync')

# Using HDF5 to handle the complex data storage without using CSV/FITS (except for the map)
store = pd.HDFStore("data/galex.h5")
store['data'] = df

# Get from the HDF5 rather than DB
df = pd.read_hdf("data/galex.h5", "data")

# Make and Store
hp_map, ptab = make_map(df)
store['ptab'] = ptab
output_map(hp_map, outfile='data/galex_map.fits')  # in FITS format
store.close()

# Read
hp_map = read_map(mapfile='data/galex_map.fits')

# Plot
make_plot(hp_map, outfile='image/mast_galex_map.png', title='GALEX')
# endregion


# region PS1
df = get_db_data(mission='PS1', limit=None,
                 constraints="AND dataproduct_type='image' AND calib_level=3 ",
                 server='mastdbtest', db='CAOMv230ProdSync')

# Using HDF5 to handle the complex data storage without using CSV/FITS (except for the map)
store = pd.HDFStore("data/ps1.h5")
store['data'] = df

# Get from the HDF5 rather than DB
# df = pd.read_hdf("data/ps1.h5", "data")

hp_map, ptab = make_map(df)

store['ptab'] = ptab
output_map(hp_map, outfile='data/ps1_map.fits')  # in FITS format

hp_map = read_map(mapfile='data/ps1_map.fits')
store.close()

make_plot(hp_map, outfile='image/mast_ps1_map.png', title='PS1')

# endregion


# Testing changes
from mast_plot import *
df = read_file_data('data/tess_data.csv', fmt='pandas')
hp_map, ptab = make_map(df)
make_plot(hp_map, outfile='image/mast_tess_map_test.png', title='TESS')

# df = read_file_data('data/hst_data.csv')
# hp_map, ptab = make_map(df)
# make_plot(hp_map, outfile='image/mast_hst_map_test.png', title='HST')

# Using HDF5 to store data
# https://pandas.pydata.org/pandas-docs/stable/user_guide/io.html#io-hdf5
store = pd.HDFStore("data/tess.h5")
store["data"] = df
store["ptab"] = ptab
store.close()

df = pd.read_hdf("data/tess.h5", "data")
ptab = pd.read_hdf("data/tess.h5", "ptab")

# Script to test things out
from db_connect import get_db_data
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
                  server='mastdbtest', db='CAOMv240ProdSync')

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
                 constraints="AND dataproduct_type='image' ",
                 server='mastdbtest', db='CAOMv240ProdSync')

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
                 server='mastdbtest', db='CAOMv240ProdSync')

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

# region Kepler
# No FFI data for Kepler yet
# TODO: create fake minimalist entry for Kepler?
# df1 = get_db_data(mission='Kepler', limit=None,
#                  constraints="AND dataproduct_type='cube' AND calib_level=2 ",
#                  server='mastdbtest', db='CAOMv230ProdSync')
df = get_db_data(mission='K2', limit=None,
                 constraints="AND dataproduct_type='image' AND calib_level=3 ",
                 server='mastdbtest', db='CAOMv240ProdSync')

# df = pd.concat([df1, df2])

# Using HDF5 to handle the complex data storage without using CSV/FITS (except for the map)
store = pd.HDFStore("data/kepler.h5")
store['data'] = df

# Get from the HDF5 rather than DB
# df = pd.read_hdf("data/ps1.h5", "data")

hp_map, ptab = make_map(df)

store['ptab'] = ptab
output_map(hp_map, outfile='data/kepler_map.fits')  # in FITS format

hp_map = read_map(mapfile='data/kepler_map.fits')
store.close()

make_plot(hp_map, outfile='image/mast_kepler_map.png', title='K2')
# endregion

# region JWST
# Ignoring planned observations and some problematic cases
df = get_db_data(mission='JWST', limit=None,
                 server='mastdbtest', db='CAOMv240ProdSync',
                 constraints="AND calib_level !=-1 " \
                    "AND proposal_id NOT IN ('1128','1678','1492','1537','1538','1231', '2344', '1132', '1125')" \
                    "AND obs_id NOT IN ('jw01189-o001_s00004_nirspec_f100lp-g140m-s200a1'," \
                    "'jw01492-o005_s00004_nirspec_f100lp-g140m-s1600a1-sub2048'," \
                    "'jw01128-o016_s00004_nirspec_f100lp-g140m-s1600a1'," \
                    "'jw01128-o004_s00004_nirspec_f100lp-g140m-s1600a1-allslits'," \
                    "'jw01536-o002_s00004_nirspec_f100lp-g140m-s1600a1-sub2048'," \
                    "'jw01538-o160_s00004_nirspec_f100lp-g140m-s1600a1-sub2048'," \
                    "'jw01189-o001_s00004_nirspec_f290lp-g395m-s200a1'," \
                    "'jw01947-o012_s00004_nirspec_f290lp-g395m-s1600a1'," \
                    "'jw01536-o002_s00004_nirspec_f170lp-g235m-s1600a1-sub2048'," \
                    "'jw01536-o002_s00004_nirspec_f290lp-g395m-s1600a1-sub2048'," \
                    "'jw01189-o001_s00004_nirspec_f170lp-g235m-s200a1'," \
                    "'jw01967-o012_s00004_nirspec_f290lp-g395m-s200a2'," \
                    "'jw01222-o002_s00004_nirspec_f170lp-g235h-s200a2'," \
                    "'jw01967-o011_s00004_nirspec_f290lp-g395m-s200a2'," \
                    "'jw01536-o002_s00004_nirspec_f070lp-g140h-s1600a1-sub2048'," \
                    "'jw01967-o005_s00004_nirspec_f290lp-g395m-s200a2'," \
                    "'jw01240-o002_s00004_nirspec_f290lp-g395m-s200a1'," \
                    "'jw01222-o002_s00004_nirspec_f070lp-g140h-s200a2'," \
                    "'jw01967-o002_s00004_nirspec_f290lp-g395m-s200a2'," \
                    "'jw01967-o003_s00004_nirspec_f290lp-g395m-s200a2'," \
                    "'jw02285-o002_s00004_nirspec_f290lp-g395m-s200a1'," \
                    "'jw01536-o002_s00004_nirspec_f170lp-g235h-s1600a1-sub2048'," \
                    "'jw01222-o012_s00004_nirspec_f170lp-g235h-s200a2'," \
                    "'jw01222-o012_s00004_nirspec_f070lp-g140h-s200a2'," \
                    "'jw01222-o002_s00004_nirspec_f070lp-g140h-s200a1'," \
                    "'jw01222-o012_s00004_nirspec_f070lp-g140h-s200a1'," \
                    "'jw01536-o002_s00004_nirspec_f290lp-g395h-s1600a1-sub2048'," \
                    "'jw01536-o002_s00004_nirspec_f100lp-g140h-s1600a1-sub2048'," \
                    "'jw01539-o076_s00004_nirspec_f100lp-g140h-s1600a1-sub2048'," \
                    "'jw01539-o003_s00004_nirspec_f100lp-g140h-s1600a1-sub2048'," \
                    "'jw01222-o012_s00004_nirspec_f170lp-g235h-s200a1'," \
                    "'jw01309-o015_s00004_nirspec_f290lp-g395h-s200a2-allslits'," \
                    "'jw01222-o002_s00004_nirspec_f170lp-g235h-s200a1'," \
                    "'jw01240-o004_s00004_nirspec_f290lp-g395h-s200a1'," \
                    "'jw02285-o002_s00004_nirspec_f170lp-g235m-s200a1'," \
                    "'jw01309-o022_s00004_nirspec_f290lp-g395h-s200a2-allslits'," \
                    "'jw01118-o002_s00004_nirspec_f100lp-g140h-s200a1-allslits'," \
                    "'jw01189-o011_s00004_nirspec_f290lp-g395h-s200a1'," \
                    "'jw01118-o002_s00004_nirspec_f100lp-g140h-s200a2-allslits'," \
                    "'jw01309-o015_s00003_nirspec_f290lp-g395h-s200a2-allslits'," \
                    "'jw01309-o022_s00003_nirspec_f290lp-g395h-s200a2-allslits'," \
                    "'jw01240-o004_s00003_nirspec_f290lp-g395h-s200a1'," \
                    "'jw01189-o011_s00003_nirspec_f290lp-g395h-s200a1'," \
                    "'jw01536-o002_s00004_nirspec_f070lp-g140m-s1600a1-sub2048'," \
                    "'jw01537-o007_s00004_nirspec_f100lp-g140m-s1600a1-sub2048')")

# Using HDF5 to handle the complex data storage without using CSV/FITS (except for the map)
store = pd.HDFStore("data/jwst.h5")
store['data'] = df

# Get from the HDF5 rather than DB
# df = pd.read_hdf("data/jwst.h5", "data")

hp_map, ptab = make_map(df)

# Identify overly large/incorrect footprints
# ptab[ptab['np'] > 1000].sort_values('np', ascending=False)
# for i, row in ptab[ptab['np'] > 100].sort_values('np', ascending=False).iterrows():
#     print(row['obs_id'])

store['ptab'] = ptab
output_map(hp_map, outfile='data/jwst_map.fits')  # in FITS format

hp_map = read_map(mapfile='data/jwst_map.fits')
store.close()

make_plot(hp_map, outfile='image/mast_jwst_map.png', title='JWST')
make_plot(hp_map, outfile='image/mast_jwst_map_v2.png', title='JWST', grids=False)
# endregion

# Testing changes
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

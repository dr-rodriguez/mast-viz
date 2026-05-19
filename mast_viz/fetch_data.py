# Module to fetch data from the database or HDF5 storage
import os
import pandas as pd
from mast_viz.utils.db_connect import get_db_data
from mast_viz.utils.mast_plot import make_plot, make_map, output_map, read_map

# --- MISSION CONFIGURATION ---
# Tweak these variables for a single mission run. 
# You can uncomment/comment these as needed.
# Use "! =" not ""= / =" for sql queries"

# MISSION = "HST"
# CONSTRAINTS = "AND dataproduct_type='image' AND provenance_name!='HAP' "

# MISSION = "TESS"
# CONSTRAINTS = "AND dataproduct_type='image' "

# MISSION = "JWST"
# CONSTRAINTS = ("AND calib_level !=-1 "
#                "AND proposal_id NOT IN ('1128','1678','1492','1537','1538','1231', '2344', '1132', '1125') ")

# MISSION = "GALEX"
# CONSTRAINTS = "AND dataproduct_type='image' "

# MISSION = "PS1"
# CONSTRAINTS = "AND dataproduct_type='image' AND calib_level=3 "

# MISSION = "K2"
# CONSTRAINTS = "AND dataproduct_type='image' AND calib_level=3 "

# MISSION = "SDSS"
# CONSTRAINTS = "AND dataproduct_type='image' "

MISSION = "HLSP"
CONSTRAINTS = ""

# Default (currently set for HST)
# MISSION = "HST"
# CONSTRAINTS = "AND dataproduct_type='image' AND provenance_name NOT IN ('HAP', 'HASP', 'HSLA') "

QUERY_FRESH = False  # Set to True to bypass HDF5 cache and query the database
RESUME = True       # Set to True to resume from temporary CSV chunks
DATA_DIR = "data"
MAKE_PLOTS = True
RUN_CHUNKS = False
NUM_CHUNKS = 20

# -----------------------------


def fetch_mission_data(mission=None, constraints=None, query_fresh=None, data_dir=None, make_plots=None, run_chunks=None, num_chunks=None, resume=None):
    """
    Fetch data for a specific mission from the database or HDF5 cache.
    
    Parameters
    ----------
    mission : str, optional
        Mission name (e.g., "HST", "JWST"). Defaults to MISSION variable.
    constraints : str, optional
        SQL constraints. Defaults to CONSTRAINTS variable.
    query_fresh : bool, optional
        Whether to bypass cache and query DB. Defaults to QUERY_FRESH variable.
    data_dir : str, optional
        Directory for data files. Defaults to DATA_DIR variable.
    make_plots : bool, optional
        Whether to generate coverage plots. Defaults to MAKE_PLOTS variable.
    run_chunks : bool, optional
        Whether to run query in chunks. Defaults to RUN_CHUNKS variable.
    num_chunks : int, optional
        Number of chunks to use. Defaults to NUM_CHUNKS variable.
    resume : bool, optional
        Whether to resume from temporary CSV chunks. Defaults to RESUME variable.
        
    Returns
    -------
    pd.DataFrame
        Fetched data.
    """
    # Use provided arguments or fall back to module-level configuration
    mission = mission or MISSION
    constraints = constraints if constraints is not None else CONSTRAINTS
    query_fresh = query_fresh if query_fresh is not None else QUERY_FRESH
    data_dir = data_dir or DATA_DIR
    make_plots = make_plots if make_plots is not None else MAKE_PLOTS
    run_chunks = run_chunks if run_chunks is not None else RUN_CHUNKS
    num_chunks = num_chunks if num_chunks is not None else NUM_CHUNKS
    resume = resume if resume is not None else RESUME
    
    h5_path = os.path.join(data_dir, f"{mission.lower()}.h5")
    fits_path = os.path.join(data_dir, f"{mission.lower()}_map.fits")
    df = None
    hp_map = None
    
    # Try to load from HDF5 if not querying fresh
    if not query_fresh and os.path.exists(h5_path):
        print(f"Loading data for {mission} from {h5_path}")
        try:
            with pd.HDFStore(h5_path, mode='r') as store:
                if 'data' in store:
                    df = store['data']
                else:
                    print(f"'data' key not found in {h5_path}. Fetching from DB.")
                
                ptab = None
                if 'ptab' in store:
                    ptab = store['ptab']
                    print(f"Loaded ptab for {mission} from cache.")
        except Exception as e:
            print(f"Error reading HDF5 for {mission}: {e}. Falling back to DB.")

    if df is not None and hp_map is None and os.path.exists(fits_path):
        try:
            hp_map = read_map(fits_path)
            print(f"Loaded map for {mission} from {fits_path}")
        except Exception as e:
            print(f"Error reading FITS map for {mission}: {e}")

    if df is None:
        # Fetch from database
        print(f"Fetching data for {mission} from database...")
        df = get_db_data(mission=mission, constraints=constraints, run_chunks=run_chunks, num_chunks=num_chunks, resume=resume)
        
        # Apply mission-specific post-processing if needed (matching script.py logic)
        if mission == "SDSS":
            bad_tmin = 40587  # 1970, this is a bug in the values
            df = df[df['t_min'] != bad_tmin]
        elif mission == "HST":
            # Remove data before 1990
            df = df[df['t_min'] >= 48005.]

        # Generate HEALPix map and pixel table (ptab)
        print(f"Generating map and ptab for {mission}...")
        hp_map, ptab = make_map(df)

        # Save to HDF5 (always update cache when fetching fresh or if cache was missing)
        print(f"Saving data and ptab for {mission} to {h5_path}")
        os.makedirs(data_dir, exist_ok=True)
        try:
            with pd.HDFStore(h5_path, mode='a') as store:
                store["data"] = df
                store["ptab"] = ptab
            
            # Save the HEALPix map to a FITS file
            print(f"Saving map for {mission} to {fits_path}")
            output_map(hp_map, outfile=fits_path)
        except Exception as e:
            print(f"Error saving data for {mission}: {e}")
    else:
        # df was loaded from cache. Check if ptab and fits map also exist.
        if ptab is None or hp_map is None:
             print(f"Map or ptab missing from cache for {mission}, generating...")
             hp_map, ptab = make_map(df)
             try:
                 with pd.HDFStore(h5_path, mode='a') as store:
                     store["ptab"] = ptab
                 output_map(hp_map, outfile=fits_path)
             except Exception as e:
                 print(f"Error updating map/ptab for {mission}: {e}")

    if make_plots:
        print(f"Generating plot for {mission}...")
        if hp_map is None:
            if os.path.exists(fits_path):
                hp_map = read_map(fits_path)
            else:
                hp_map, _ = make_map(df)
        
        plot_path = os.path.join("image", f"mast_{mission.lower()}_map.png")
        os.makedirs("image", exist_ok=True)
        make_plot(hp_map, outfile=plot_path, title=mission)
            
    return df


if __name__ == "__main__":
    # When run as a script, use the configuration at the top
    data = fetch_mission_data()
    print(f"Fetched {len(data)} rows for {MISSION}")
    if not data.empty:
        print(data.head())

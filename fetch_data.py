# Module to fetch data from the database or HDF5 storage
import os
import pandas as pd
from db_connect import get_db_data
from mast_plot import make_plot, make_map

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

# Default (currently set for HST)
MISSION = "HST"
CONSTRAINTS = "AND dataproduct_type='image' AND provenance_name NOT IN ('HAP', 'HASP', 'HSLA') "

QUERY_FRESH = True  # Set to True to bypass HDF5 cache and query the database
DATA_DIR = "data"
MAKE_PLOTS = False
# -----------------------------


def fetch_mission_data(mission=None, constraints=None, query_fresh=None, data_dir=None, make_plots=None):
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
    
    h5_path = os.path.join(data_dir, f"{mission.lower()}.h5")
    df = None
    
    # Try to load from HDF5 if not querying fresh
    if not query_fresh and os.path.exists(h5_path):
        print(f"Loading data for {mission} from {h5_path}")
        try:
            with pd.HDFStore(h5_path, mode='r') as store:
                if 'data' in store:
                    df = store['data']
                else:
                    print(f"'data' key not found in {h5_path}. Fetching from DB.")
        except Exception as e:
            print(f"Error reading HDF5 for {mission}: {e}. Falling back to DB.")

    if df is None:
        # Fetch from database
        print(f"Fetching data for {mission} from database...")
        df = get_db_data(mission=mission, constraints=constraints)
        
        # Apply mission-specific post-processing if needed (matching script.py logic)
        if mission == "SDSS":
            bad_tmin = 40587  # 1970, this is a bug in the values
            df = df[df['t_min'] != bad_tmin]
        elif mission == "HST":
            # Remove data before 1990
            df = df[df['t_min'] >= 48005.]

        # Save to HDF5 (always update cache when fetching fresh or if cache was missing)
        print(f"Saving data for {mission} to {h5_path}")
        os.makedirs(data_dir, exist_ok=True)
        try:
            with pd.HDFStore(h5_path, mode='a') as store:
                store["data"] = df
        except Exception as e:
            print(f"Error saving to HDF5: {e}")

    if make_plots:
        print(f"Generating plot for {mission}...")
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

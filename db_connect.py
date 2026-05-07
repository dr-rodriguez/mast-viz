# Script handling database connections
import pyodbc
import pandas as pd
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv(".env", override=True)
MAST_DB_SERVER = os.getenv("MAST_DB_SERVER")
MAST_DB_NAME = os.getenv("MAST_DB_NAME")


def get_db_data(mission="HST", constraints="", limit=None):
    # Get the data from ObsPointing, return as pandas dataframe
    # Fields would be observationID, exposure time, wavelength range, s_region

    if limit:
        limit_string = f"TOP {limit}"
    else:
        limit_string = ""

    sql = f"""SELECT {limit_string} 
        obs_id, obs_collection, instrument_name, dataproduct_type, s_region, t_exptime, t_min, em_min, em_max
        from ObsPointing with (nolock) 
        WHERE obs_collection='{mission}' AND t_exptime IS NOT NULL AND t_min IS NOT NULL AND s_region IS NOT NULL 
        {constraints} 
        ORDER BY t_min"""

    conn = pyodbc.connect(dsn=MAST_DB_SERVER, database=MAST_DB_NAME, trusted_connection="yes", autocommit=True)
    df = pd.read_sql(sql, conn)

    return df

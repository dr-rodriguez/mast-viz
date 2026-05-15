# Script handling database connections
import pyodbc
import pandas as pd
import os
from dotenv import load_dotenv

# Load environment variables
load_dotenv(".env", override=True)
MAST_DB_SERVER = os.getenv("MAST_DB_SERVER")
MAST_DB_NAME = os.getenv("MAST_DB_NAME")


def get_db_data(mission="HST", constraints="", limit=None, run_chunks=False, num_chunks=10):
    # Get the data from ObsPointing, return as pandas dataframe
    # Fields would be observationID, exposure time, wavelength range, s_region

    def build_sql(limit_val=None, extra_constraints=""):
        limit_string = f"TOP {limit_val}" if limit_val else ""
        return f"""SELECT {limit_string} 
            obs_id, obs_collection, instrument_name, dataproduct_type, s_region, t_exptime, t_min, em_min, em_max
            from ObsPointing with (nolock) 
            WHERE obs_collection='{mission}' AND t_exptime IS NOT NULL AND t_min IS NOT NULL AND s_region IS NOT NULL 
            {constraints} {extra_constraints}"""

    print("Connecting to DB...")
    conn = pyodbc.connect(dsn=MAST_DB_SERVER, database=MAST_DB_NAME, trusted_connection="yes", autocommit=True)

    if run_chunks:
        print(f"Running in {num_chunks} chunks...")
        range_sql = f"""SELECT MIN(obsid) as min_obsid, MAX(obsid) as max_obsid 
        FROM ObsPointing WITH (NOLOCK) 
        WHERE obs_collection='{mission}' 
        AND t_exptime IS NOT NULL AND t_min IS NOT NULL AND s_region IS NOT NULL 
        {constraints}"""
        range_df = pd.read_sql(range_sql, conn)
        min_id = range_df["min_obsid"].iloc[0]
        max_id = range_df["max_obsid"].iloc[0]

        if min_id is None or max_id is None:
            print("No data found for range.")
            return pd.DataFrame()

        step = (max_id - min_id) // num_chunks + 1
        dfs = []
        for i in range(num_chunks):
            start = min_id + i * step
            end = start + step - 1
            if i == num_chunks - 1:
                end = max_id
            
            print(f"  Chunk {i+1}/{num_chunks}: obsid {start} to {end}")
            chunk_sql = build_sql(limit, f"AND obsid BETWEEN {start} AND {end}")
            dfs.append(pd.read_sql(chunk_sql, conn))
        
        df = pd.concat(dfs, ignore_index=True)
    else:
        print("Running single query...")
        sql = build_sql(limit)
        df = pd.read_sql(sql, conn)

    # Sort by t_min
    print("Sorting by t_min...")
    df = df.sort_values(by="t_min", ascending=True).reset_index()

    print(f"Returning {len(df)} rows.")

    return df

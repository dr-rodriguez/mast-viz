import matplotlib

matplotlib.use("Agg")
import pandas as pd
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
import os
from astropy.table import Table
from matplotlib import cm
import copy
import healpy.rotator as R
from mast_viz.utils.utils import parse_s_region, get_polygon_cdshealpix

plt.interactive(False)


def save_db_data(df, datafile="data.csv", fmt="fits"):
    if fmt.lower() == "csv":
        # Output as CSV file
        df.to_csv(datafile, index=False)
    else:
        # Output as FITS file
        t2 = Table.from_pandas(df, index=False)
        t2.write(datafile, overwrite=True)


def read_file_data(datafile="data.fits", fmt="table"):
    # Read the data, either as a csv into pandas or fits into astropy Table
    if datafile.endswith("fits"):
        t = Table.read(datafile)
    else:
        t = Table.read(datafile, format="csv")
        # df = pd.read_csv(datafile)

    # Decode bytes columns
    for col in t.colnames:
        if isinstance(t[col][0], bytes):
            t[col] = [x.decode() for x in t[col]]

    if fmt.lower() == "table":
        return t
    else:
        return t.to_pandas()


def _resume_paths(mission=None):
    """
    Return HEALPix resume file paths and a log label.

    When mission is omitted, use shared cross-mission temp files.
    """
    os.makedirs("data", exist_ok=True)
    if mission is None:
        return "data/temp_hpmap.csv", "data/temp_ptab.csv", "cross-mission"
    tag = mission.lower()
    return f"data/temp_{tag}_hpmap.csv", f"data/temp_{tag}_ptab.csv", mission


def make_map(df, nside=256, exp_col="t_exptime", verbose=False, mission=None):
    """
    Make Healpix map

    Parameters
    ----------
    df
        Source dataframe
    nside
        HEALPix resolution
    exp_col
        Exposure column name
    mission
        Mission name for mission-specific resume files. When omitted, resume
        state is stored in shared ``data/temp_hpmap.csv`` and
        ``data/temp_ptab.csv`` for cross-mission maps.

    Returns
    -------
    hp_map
        Exposure map by HEALPix
    ptab
        HEALPix pixel values for the provided dataframe
    """
    temp_hpmap_path, temp_ptab_path, resume_label = _resume_paths(mission=mission)

    # number of pixels for that resolution
    npix = hp.nside2npix(nside)
    resolution = hp.nside2resol(nside, arcmin=True)
    print(f"NSIDE={nside} NPIX={npix} Resolution(arcmin)={resolution}")

    # Try to resume from temporary files
    hp_map = None
    ptab = []
    processed_indices = set()

    if os.path.exists(temp_hpmap_path) and os.path.exists(temp_ptab_path):
        print(f"Temporary files found for {resume_label}. Attempting to resume...")
        try:
            hp_map_df = pd.read_csv(temp_hpmap_path)
            if len(hp_map_df) == npix:
                hp_map = hp_map_df["value"].values.copy()
                
                ptab_df = pd.read_csv(temp_ptab_path)
                if not ptab_df.empty:
                    processed_indices = set(ptab_df["i"].values)
                    ptab = ptab_df.to_dict("records")
                    # Fix 'ind' column which is stored as string in CSV
                    for d in ptab:
                        if isinstance(d["ind"], str):
                            # Remove brackets and split by space/newline
                            s = d["ind"].strip("[]").replace("\n", " ")
                            d["ind"] = np.fromstring(s, sep=" ", dtype=int)
                    print(f"Resuming from row {len(processed_indices)}")
                else:
                    hp_map = np.zeros(npix)
                    ptab = []
            else:
                print(f"Temp map size {len(hp_map_df)} does not match npix {npix}. Starting fresh.")
        except Exception as e:
            print(f"Error loading temporary files: {e}. Starting fresh.")

    if hp_map is None:
        hp_map = np.zeros(npix)
        ptab = []
        processed_indices = set()

    # I believe this is an estimate for how many HP-resolution elements there are in an area,
    # which for us is very variable. I'm using the TESS FFI size of just over 12 degrees
    # However, this will lead to gigantic fits as the way the output is produced is to store
    # zeros for all possible ranges
    # We may want to construct a sparser data structure which we can then read and interpret
    size_deg = 13  # largest area to consider, in degrees
    nind = int(size_deg**2 * np.pi / (resolution / 60) ** 2 * size_deg)

    # ptab = np.zeros(len(df), dtype=[('np', 'i8'), ('ind', '%di8' % nind)])  # table to store healpix ids
     
    length = len(df)
    count = 0

    for i, row in df.iterrows():
        if i in processed_indices:
            continue

        # Print status every 10%
        if count % (length // 10 + 1) == 0:
            print(f"Progress: {count}/{length} ({count/length*100:.1f}%)")
        
        count += 1
        
        if row["s_region"] is not None and row[exp_col] is not None:
            # Parse the footprint
            try:
                coords = parse_s_region(row["s_region"])  # converts CIRCLE to 16-point POLYGON
                ra_list, dec_list = coords["ra"], coords["dec"]
            except:
                print("Unable to parse s_region for {}".format(row["obs_id"]))
                continue

            # Generate HEALPix indices for the footprint
            try:
                # CDSHEALPIX package
                # Indices of all pixels inside or intersecting the polygon
                ipix = get_polygon_cdshealpix(ra_list, dec_list, depth=int(np.log2(nside)))
                if verbose:
                    print("CDS", row["obs_id"], ipix)
            except:
                # HEALPY package
                try:
                    # Indices of all pixels inside or intersecting the convex polygon vec (if not convex it will fail)
                    vec = hp.ang2vec(ra_list, dec_list, lonlat=True)
                    ipix = hp.query_polygon(nside, vec, inclusive=True)
                    if verbose:
                        print("HP ", row["obs_id"], ipix)
                except:
                    print("Unable to get HP indeces for {}: {}".format(row["obs_id"], row["s_region"]))
                    continue

            # Populate output lists
            try:
                hp_map[ipix] = hp_map[ipix] + float(row[exp_col])  # adding up exposures for that healpix pixel

                p_dict = {
                    "i": i,
                    "obs_id": row["obs_id"],  # identifiers
                    "np": len(ipix),  # number of pixels for each observation
                    "ind": ipix,
                }  # pixel values for each observation
                ptab.append(p_dict)
                # ptab['np'][i] = len(ipix)  # number of pixels for each observation
                # ptab['ind'][i][0:len(ipix)] = ipix  # pixel values for each observation
            except:
                print("Unable to store indices for {}: {}".format(row["obs_id"], ipix))
                continue

        # Periodically save state (every 1000 items)
        if count % 1000 == 0:
            pd.DataFrame({"value": hp_map}).to_csv(temp_hpmap_path, index=False)
            pd.DataFrame(ptab).to_csv(temp_ptab_path, index=False)

    # Save final state before returning
    pd.DataFrame({"value": hp_map}).to_csv(temp_hpmap_path, index=False)
    pd.DataFrame(ptab).to_csv(temp_ptab_path, index=False)

    pdf = pd.DataFrame(ptab)

    return hp_map, pdf


def output_map(hp_map, outfile="mast_map.fits"):
    # Output the map
    hp.write_map(outfile, hp_map, coord="C", overwrite=True)
    return


def read_map(mapfile):
    return hp.read_map(mapfile)


def _project_lonlat(ax, lon, lat, coord="C"):
    """Project lon/lat to axes coordinates if the point lies in the plot area."""
    vec = R.dir2vec(lon, lat, lonlat=True)
    vec = R.Rotator(coord=ax.proj.mkcoord(coord=coord)[::-1]).I(vec)
    x, y = ax.proj.vec2xy(vec, direct=False)
    if not np.isfinite(x) or not np.isfinite(y):
        return None

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    pad = 0.03 * max(xmax - xmin, ymax - ymin)
    if x < xmin + pad or x > xmax - pad or y < ymin + pad or y > ymax - pad:
        return None
    return float(x), float(y)


def _add_graticule_labels(ax, dpar=45.0, dmer=30.0, coord="C", fontsize=8, color="white"):
    """Add small RA/Dec labels on visible major graticule lines."""
    for ra in np.arange(0.0, 360.0, dmer):
        if _project_lonlat(ax, ra, 0.0, coord=coord) is not None:
            ax.projtext(
                ra,
                0.0,
                f"{int(ra)}°",
                lonlat=True,
                coord=coord,
                color=color,
                fontsize=fontsize,
                ha="center",
                va="bottom",
            )

    for dec in np.arange(-90.0 + dpar, 90.0, dpar):
        if abs(dec) < 1e-9:
            continue

        best_ra = None
        best_score = None
        for ra in np.linspace(0.0, 360.0, 144, endpoint=False):
            xy = _project_lonlat(ax, ra, dec, coord=coord)
            if xy is None:
                continue
            score = abs(xy[0])
            if best_score is None or score > best_score:
                best_score = score
                best_ra = ra

        if best_ra is not None:
            dec_label = f"+{int(dec)}°" if dec > 0 else f"{int(dec)}°"
            ax.projtext(
                best_ra,
                dec,
                dec_label,
                lonlat=True,
                coord=coord,
                color=color,
                fontsize=fontsize,
                ha="center",
                va="center",
            )


def make_plot(
    hp_map,
    outfile="mast_map.png",
    title="",
    dpi=300,
    grids=True,
    grid_labels=False,
    grid_label_fontsize=8,
):
    # Generate the map

    SKYCOLOR = '#003B4D'  # MAST darkest turquoise

    pngfile2 = os.path.splitext(outfile)[0] + f"_{dpi}.png"

    # Plot options
    plt.style.use("dark_background")
    cmap = copy.copy(cm.get_cmap("cividis"))
    cmap.set_bad(SKYCOLOR)
    cmap.set_under('k')

    plt.rcParams.update({"font.size": 15})
    lon = np.arange(360)
    lat = np.zeros(360)

    # datestr = time.strftime('%B %d, %Y')
    hp.mollview(
        np.log10(hp_map + 0.1),
        cmap=cmap,
        min=0.0,
        max=6.0,  # exptime limits
        rot=-80,
        flip="geo",
        coord="C",
        cbar=False,  # color bar
        notext=True,
        title=title,
        bgcolor="black",
        badcolor=SKYCOLOR,
        norm="linear",
        xsize=1000,
    )

    if grids:
        hp.projplot(lon, lat, "r", lonlat=True, coord="G")
        hp.graticule(dpar=45.0, dmer=30.0, coord="C", color="white")
        if grid_labels:
            _add_graticule_labels(
                plt.gca(),
                dpar=45.0,
                dmer=30.0,
                coord="C",
                fontsize=grid_label_fontsize,
                color="white",
            )

    plt.savefig(pngfile2, dpi=dpi)

    plt.close()

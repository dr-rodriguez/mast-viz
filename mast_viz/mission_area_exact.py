"""
Compute sky coverage from exact s_region footprints in cached mission HDF5 files.

Reads data/<mission>.h5 (/data key) and unions STC-S POLYGON/CIRCLE footprints into
a MOC at a per-mission HEALPix order. More accurate than mission_area.py, which
quantizes to nside=256 and counts intersecting pixels.

Per-mission default orders (override with --order):
  hst/jwst: 18, sdss/galex/hlsp: 16, ps1/kepler: 14, tess: 11

HEALPix cell area at each order (nside = 2^order; see ORDER_CELL_DEG2 below):
  order 11 (tess):      ~8.20e-4 deg^2 (~2.9 arcmin)
  order 14 (ps1, ...):  ~1.28e-5 deg^2 (~0.36 arcmin)
  order 16 (sdss, ...): ~8.00e-7 deg^2 (~0.09 arcmin)
  order 18 (hst, jwst): ~5.00e-8 deg^2 (~0.02 arcmin)
  order 8 (mission_area.py nside=256): ~0.0525 deg^2 (~11.5 arcmin)

TESS and PS1 use shallow orders because their footprints span degrees; deep orders
are prohibitively slow and unnecessary for boundary accuracy on large tiles.

Notes
-----
- Full sky reference: 4*pi sr = 41252.96 deg^2.
- Cached MOCs are written to data/<mission>_moc_o<order>.fits.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from collections import defaultdict

import astropy.units as u
import healpy as hp
import numpy as np
import pandas as pd
from mocpy import MOC

# Full-sky area constants used to convert MOC sky_fraction to deg^2 and sr.
FULL_SKY_DEG2 = 4.0 * np.pi * (180.0 / np.pi) ** 2
FULL_SKY_SR = 4.0 * np.pi

DEFAULT_ORDER = 16
# Per-mission HEALPix order defaults: tuned to footprint angular scale.
MISSION_ORDER: dict[str, int] = {
    "hst": 18,
    "jwst": 18,
    "sdss": 16,
    "galex": 16,
    "hlsp": 16,
    "ps1": 14,
    "kepler": 14,
    "tess": 11,
}

# HEALPix cell area (deg^2) at each order used above (nside = 2^order).
# MOC boundary error is roughly one cell width along each footprint edge.
# mission_area.py uses nside=256 (order 8), ~1000x coarser than hst/jwst order 18.
MAP_BASED_ORDER = 8  # nside=256 in mission_area.py / make_map
_ORDER_REFERENCE_ORDERS = frozenset(MISSION_ORDER.values()) | {DEFAULT_ORDER, MAP_BASED_ORDER}
ORDER_CELL_DEG2: dict[int, float] = {
    order: hp.nside2pixarea(2 ** order, degrees=True)
    for order in _ORDER_REFERENCE_ORDERS
}


def available_missions(data_dir: str = "data") -> list[str]:
    """Return mission names with cached HDF5 files, sorted."""
    pattern = os.path.join(data_dir, "*.h5")
    missions = []
    for path in glob.glob(pattern):
        basename = os.path.basename(path)
        if basename.endswith(".h5"):
            missions.append(basename[:-len(".h5")])
    return sorted(missions)


def default_order(mission: str) -> int:
    """Return the default HEALPix order for a mission."""
    return MISSION_ORDER.get(mission.lower(), DEFAULT_ORDER)


def moc_cache_path(data_dir: str, mission: str, order: int) -> str:
    """Path for a cached mission MOC FITS file."""
    return os.path.join(data_dir, f"{mission.lower()}_moc_o{order}.fits")


def is_cache_valid(cache_path: str, h5_path: str) -> bool:
    """True when the cache exists and is newer than the source HDF5."""
    return (
        os.path.exists(cache_path)
        and os.path.exists(h5_path)
        and os.path.getmtime(cache_path) >= os.path.getmtime(h5_path)
    )


def parse_regions(
    regions: np.ndarray,
) -> tuple[dict[int, list[tuple[np.ndarray, np.ndarray]]], tuple | None, int]:
    """
    Parse s_region STC-S strings into polygon groups and circle arrays.

    Returns
    -------
    poly_groups
        Dict keyed by vertex count; values are (lon, lat) array pairs.
    cones
        (lon, lat, radius) float arrays, or None when no circles.
    n_unsupported
        Count of rows that could not be parsed.
    """
    # Polygons grouped by vertex count (mocpy requires uniform polygons per batch).
    poly_groups: dict[int, list[tuple[np.ndarray, np.ndarray]]] = defaultdict(list)
    cone_lons: list[float] = []
    cone_lats: list[float] = []
    cone_radii: list[float] = []
    n_unsupported = 0

    for reg in regions:
        if reg is None:
            n_unsupported += 1
            continue

        if isinstance(reg, bytes):
            reg = reg.decode()

        # Split STC-S into tokens; first token is shape type (POLYGON/CIRCLE).
        toks = str(reg).strip().split()
        if not toks:
            n_unsupported += 1
            continue

        head = toks[0].upper()
        # Extract numeric coordinates only; skip labels like ICRS.
        nums: list[float] = []
        for token in toks[1:]:
            try:
                nums.append(float(token))
            except ValueError:
                continue

        if head == "POLYGON":
            if len(nums) < 6:
                n_unsupported += 1
                continue
            arr = np.asarray(nums)
            # Interleaved lon/lat pairs; normalize lon to [0, 360) for mocpy.
            lon = arr[0::2] % 360.0
            lat = arr[1::2]
            if len(lon) != len(lat) or len(lon) < 3:
                n_unsupported += 1
                continue
            # Drop repeated closing vertex when first == last.
            if (
                abs(lon[0] - lon[-1]) < 1e-9
                and abs(lat[0] - lat[-1]) < 1e-9
            ):
                lon = lon[:-1]
                lat = lat[:-1]
            if len(lon) < 3:
                n_unsupported += 1
                continue
            poly_groups[len(lon)].append((lon, lat))
        elif head == "CIRCLE":
            if len(nums) < 3:
                n_unsupported += 1
                continue
            # Center lon/lat and radius in degrees (exact circles, not polygons).
            cone_lons.append(nums[0] % 360.0)
            cone_lats.append(nums[1])
            cone_radii.append(nums[2])
        else:
            # RANGE, UNION, malformed rows, etc.
            n_unsupported += 1

    cones = None
    if cone_lons:
        cones = (
            np.asarray(cone_lons, dtype=np.float64),
            np.asarray(cone_lats, dtype=np.float64),
            np.asarray(cone_radii, dtype=np.float64),
        )

    return poly_groups, cones, n_unsupported


def _union_batch(mocs: list[MOC]) -> MOC | None:
    """Union a list of MOCs; return None for an empty list."""
    if not mocs:
        return None
    if len(mocs) == 1:
        return mocs[0]
    return mocs[0].union(*mocs[1:])


def _union_chunked(mocs: list[MOC], chunk_size: int) -> MOC | None:
    """Union many MOCs in chunks to limit peak memory."""
    if not mocs:
        return None

    # Union footprints in batches, then union the batch results.
    chunk_results: list[MOC] = []
    for start in range(0, len(mocs), chunk_size):
        batch = mocs[start : start + chunk_size]
        chunk_moc = _union_batch(batch)
        if chunk_moc is not None:
            chunk_results.append(chunk_moc)

    return _union_batch(chunk_results)


def build_moc(
    poly_groups: dict[int, list[tuple[np.ndarray, np.ndarray]]],
    cones: tuple | None,
    order: int,
    chunk_size: int = 50_000,
    progress: bool = True,
) -> MOC:
    """
    Build a union MOC from parsed polygon groups and optional circle arrays.
    """
    running: MOC | None = None
    pending: list[MOC] = []
    n_poly = sum(len(items) for items in poly_groups.values())

    if progress and n_poly:
        print(
            f"Building MOC from {n_poly:,} polygons at order {order}...",
            file=sys.stderr,
        )

    # Step 1: convert each polygon group to individual footprint MOCs.
    for n_vertices, items in sorted(poly_groups.items()):
        flat: list[np.ndarray] = []
        for lon, lat in items:
            flat.append(lon)
            flat.append(lat)
        # One MOC per polygon; all polygons in this group share vertex count.
        group_mocs = MOC.from_polygons(
            np.array(flat, dtype=np.float64),
            max_depth=order,
        )
        pending.extend(group_mocs)

        # Step 2: periodically union pending footprints into the running MOC.
        if len(pending) >= chunk_size:
            chunk_moc = _union_chunked(pending, chunk_size)
            if chunk_moc is not None:
                running = (
                    chunk_moc if running is None else running.union(chunk_moc)
                )
            pending = []

    # Union any remaining polygon MOCs not yet folded into running.
    if pending:
        chunk_moc = _union_chunked(pending, chunk_size)
        if chunk_moc is not None:
            running = chunk_moc if running is None else running.union(chunk_moc)

    # Step 3: add circle footprints (e.g. GALEX) as exact cones, not polygons.
    if cones is not None:
        lon, lat, radius = cones
        if progress:
            print(f"Adding {len(lon):,} circles...", file=sys.stderr)
        cone_moc = MOC.from_cones(
            lon=lon * u.deg,
            lat=lat * u.deg,
            radius=radius * u.deg,
            max_depth=order,
            union_strategy="small_cones",
        )
        running = cone_moc if running is None else running.union(cone_moc)

    if running is None:
        return MOC(max_depth=order)

    return running


def _stats_from_moc(moc: MOC, order: int) -> tuple[int, float, float, float]:
    """Return n_cells, area_deg2, area_sr, sky_fraction for a MOC."""
    sky_fraction = float(moc.sky_fraction)
    area_sr = sky_fraction * FULL_SKY_SR
    area_deg2 = sky_fraction * FULL_SKY_DEG2
    n_cells = int(moc.n_cells(order))
    return n_cells, area_deg2, area_sr, sky_fraction


def mission_area_exact(
    mission: str,
    data_dir: str = "data",
    order: int | None = None,
    min_exptime: float = 0.0,
    use_cache: bool = True,
    recompute: bool = False,
    chunk_size: int = 50_000,
    progress: bool = True,
) -> dict:
    """
    Compute sky coverage from exact footprints in a mission HDF5 file.

    Parameters
    ----------
    mission : str
        Mission name (e.g. "HST", "jwst").
    data_dir : str
        Directory containing <mission>.h5 files.
    order : int, optional
        HEALPix order for the MOC. Defaults to the per-mission table.
    min_exptime : float
        Minimum t_exptime for an observation to be included.
    use_cache : bool
        Read/write cached MOC FITS when min_exptime is zero.
    recompute : bool
        Force rebuild even when a valid cache exists.
    chunk_size : int
        Chunk size for unioning footprint MOCs.
    progress : bool
        Print progress messages to stderr.

    Returns
    -------
    dict
        Coverage statistics for the mission.
    """
    mission_key = mission.lower()
    h5_path = os.path.join(data_dir, f"{mission_key}.h5")

    if not os.path.exists(h5_path):
        available = available_missions(data_dir)
        available_msg = ", ".join(available) if available else "(none)"
        raise FileNotFoundError(
            f"No HDF5 file found for mission '{mission}' at {h5_path}. "
            f"Available missions: {available_msg}"
        )

    if order is None:
        order = default_order(mission_key)

    cache_path = moc_cache_path(data_dir, mission_key, order)
    # Caching only applies to the full dataset (min_exptime == 0).
    can_cache = use_cache and min_exptime <= 0.0
    cached = False

    # Fast path: load precomputed MOC if cache is valid and not forced to rebuild.
    if can_cache and not recompute and is_cache_valid(cache_path, h5_path):
        if progress:
            print(f"Loading cached MOC from {cache_path}...", file=sys.stderr)
        moc = MOC.load(cache_path, format="fits")
        cached = True
        n_obs = 0
        n_used = 0
        n_unsupported = 0
    else:
        # Slow path: read observations and build MOC from s_region footprints.
        if progress:
            print(f"Reading {h5_path}...", file=sys.stderr)
        df = pd.read_hdf(h5_path, "data")

        if min_exptime > 0:
            df = df[df["t_exptime"] > min_exptime]
            if progress:
                print(
                    f"Filtered to {len(df):,} rows with t_exptime > {min_exptime}",
                    file=sys.stderr,
                )

        n_obs = len(df)
        poly_groups, cones, n_unsupported = parse_regions(
            df["s_region"].to_numpy()
        )
        n_used = n_obs - n_unsupported

        moc = build_moc(
            poly_groups,
            cones,
            order,
            chunk_size=chunk_size,
            progress=progress,
        )

        if can_cache:
            if progress:
                print(f"Saving MOC to {cache_path}...", file=sys.stderr)
            os.makedirs(data_dir, exist_ok=True)
            moc.save(cache_path, format="fits", overwrite=True)
            cached = True

    # Derive area from MOC sky_fraction (unique sky covered, overlaps removed).
    cell_deg2 = hp.nside2pixarea(2 ** order, degrees=True)
    n_cells, area_deg2, area_sr, sky_fraction = _stats_from_moc(moc, order)

    return {
        "mission": mission.upper(),
        "path": h5_path,
        "moc_path": cache_path if cached else None,
        "order": order,
        "cell_deg2": cell_deg2,
        "n_obs": n_obs,
        "n_used": n_used,
        "n_unsupported": n_unsupported,
        "n_cells": n_cells,
        "area_deg2": area_deg2,
        "area_sr": area_sr,
        "sky_fraction": sky_fraction,
        "cached": cached,
        "min_exptime": min_exptime,
    }


def format_result(
    result: dict,
    show_cells: bool = False,
    compare_result: dict | None = None,
) -> str:
    """Format exact coverage statistics for a single mission."""
    source = result["path"]
    if result["cached"] and result.get("moc_path"):
        source = f"{result['moc_path']} (cached)"

    lines = [
        (
            f"{result['mission']} exact coverage "
            f"({source}, order={result['order']})"
        ),
    ]

    if result["n_obs"] > 0:
        lines.append(
            f"  observations   : {result['n_used']:,} used / "
            f"{result['n_obs']:,} total"
        )
    if result["n_unsupported"]:
        lines.append(f"  unsupported    : {result['n_unsupported']:,}")
    if show_cells:
        lines.append(f"  MOC cells      : {result['n_cells']:,}")

    lines.extend([
        f"  area           : {result['area_deg2']:.2f} deg^2",
        f"                 : {result['area_sr']:.4f} sr",
        (
            f"  sky fraction   : {result['sky_fraction']:.5f} "
            f"({result['sky_fraction'] * 100:.3f}%)"
        ),
    ])

    if result["min_exptime"] > 0:
        lines.append(f"  min exptime    : {result['min_exptime']}")

    if compare_result is not None:
        map_area = compare_result["area_deg2"]
        diff = result["area_deg2"] - map_area
        pct = (diff / map_area * 100) if map_area else 0.0
        lines.append(f"  map-based area : {map_area:.2f} deg^2")
        lines.append(f"  difference     : {diff:+.2f} deg^2 ({pct:+.2f}%)")

    return "\n".join(lines)


def format_all_results(results: list[dict], show_cells: bool = False) -> str:
    """Format exact coverage statistics as a table for multiple missions."""
    if show_cells:
        header = (
            f"{'Mission':<8} {'Order':>5} {'Cells':>14} "
            f"{'deg^2':>12} {'sr':>10} {'Fraction':>10} {'Percent':>9}"
        )
    else:
        header = (
            f"{'Mission':<8} {'Order':>5} {'deg^2':>12} "
            f"{'sr':>10} {'Fraction':>10} {'Percent':>9}"
        )
    separator = "-" * len(header)
    rows = [header, separator]

    for result in results:
        if show_cells:
            rows.append(
                f"{result['mission']:<8} "
                f"{result['order']:>5} "
                f"{result['n_cells']:>14,} "
                f"{result['area_deg2']:>12.2f} "
                f"{result['area_sr']:>10.4f} "
                f"{result['sky_fraction']:>10.5f} "
                f"{result['sky_fraction'] * 100:>8.3f}%"
            )
        else:
            rows.append(
                f"{result['mission']:<8} "
                f"{result['order']:>5} "
                f"{result['area_deg2']:>12.2f} "
                f"{result['area_sr']:>10.4f} "
                f"{result['sky_fraction']:>10.5f} "
                f"{result['sky_fraction'] * 100:>8.3f}%"
            )

    return "\n".join(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Report sky coverage from exact s_region footprints in cached HDF5 files."
        )
    )
    parser.add_argument(
        "mission",
        nargs="?",
        help="Mission name (e.g. hst, jwst, tess)",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Report coverage for every mission with a cached HDF5 file",
    )
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Directory containing <mission>.h5 files (default: data)",
    )
    parser.add_argument(
        "--min-exptime",
        type=float,
        default=0.0,
        help="Minimum t_exptime for an observation to be included",
    )
    parser.add_argument(
        "--order",
        type=int,
        default=None,
        help="HEALPix order for the MOC (default: per-mission table)",
    )
    parser.add_argument(
        "--recompute",
        action="store_true",
        help="Rebuild the MOC even when a valid cache exists",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Do not read or write cached MOC files",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=50_000,
        help="Chunk size for unioning footprint MOCs (default: 50000)",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="Also report map-based coverage from mission_area.py",
    )
    parser.add_argument(
        "--cells",
        action="store_true",
        help="Include MOC cell counts in the output",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.all and args.mission:
        parser.error("Specify either a mission name or --all, not both.")
    if not args.all and not args.mission:
        parser.error("Specify a mission name or use --all.")

    use_cache = not args.no_cache

    if args.all:
        missions = available_missions(args.data_dir)
        if not missions:
            parser.error(f"No mission HDF5 files found in {args.data_dir}")
        results = [
            mission_area_exact(
                mission,
                data_dir=args.data_dir,
                order=args.order,
                min_exptime=args.min_exptime,
                use_cache=use_cache,
                recompute=args.recompute,
                chunk_size=args.chunk_size,
            )
            for mission in missions
        ]
        print(format_all_results(results, show_cells=args.cells))
        return 0

    # Single mission: compute exact area, optionally compare to map-based estimate.
    result = mission_area_exact(
        args.mission,
        data_dir=args.data_dir,
        order=args.order,
        min_exptime=args.min_exptime,
        use_cache=use_cache,
        recompute=args.recompute,
        chunk_size=args.chunk_size,
    )

    compare_result = None
    if args.compare:
        from mast_viz.mission_area import mission_area

        try:
            compare_result = mission_area(
                args.mission,
                data_dir=args.data_dir,
                min_exptime=args.min_exptime,
            )
        except FileNotFoundError as exc:
            print(f"Warning: {exc}", file=sys.stderr)

    print(format_result(result, show_cells=args.cells, compare_result=compare_result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Compute sky coverage for a MAST mission from cached HEALPix exposure maps.

Coverage is derived from data/<mission>_map.fits files produced by fetch_data.py.
Each map is a healpy RING map at nside=256 holding summed t_exptime per pixel.
A pixel is considered covered when its exposure exceeds --min-exptime.

Notes
-----
- Area is quantized to the nside=256 pixel size (~0.0525 deg^2) and is a slight
  overestimate because make_map collects pixels that merely intersect each footprint.
- Full sky reference: 4*pi sr = 41252.96 deg^2.
"""

from __future__ import annotations

import argparse
import glob
import os

import healpy as hp
import numpy as np

FULL_SKY_DEG2 = 4.0 * np.pi * (180.0 / np.pi) ** 2


def available_missions(data_dir: str = "data") -> list[str]:
    """Return mission names with cached map FITS files, sorted."""
    pattern = os.path.join(data_dir, "*_map.fits")
    missions = []
    for path in glob.glob(pattern):
        basename = os.path.basename(path)
        if basename.endswith("_map.fits"):
            missions.append(basename[: -len("_map.fits")])
    return sorted(missions)


def mission_area(
    mission: str,
    data_dir: str = "data",
    min_exptime: float = 0.0,
) -> dict:
    """
    Compute sky coverage for a mission from its cached exposure map.

    Parameters
    ----------
    mission : str
        Mission name (e.g. "HST", "jwst").
    data_dir : str
        Directory containing <mission>_map.fits files.
    min_exptime : float
        Minimum exposure time for a pixel to count as covered.

    Returns
    -------
    dict
        Coverage statistics for the mission.
    """
    mission_key = mission.lower()
    map_path = os.path.join(data_dir, f"{mission_key}_map.fits")

    if not os.path.exists(map_path):
        available = available_missions(data_dir)
        available_msg = ", ".join(available) if available else "(none)"
        raise FileNotFoundError(
            f"No map found for mission '{mission}' at {map_path}. "
            f"Available missions: {available_msg}"
        )

    hp_map = np.nan_to_num(hp.read_map(map_path))
    npix = len(hp_map)
    nside = hp.npix2nside(npix)
    n_covered = int(np.count_nonzero(hp_map > min_exptime))

    pix_sr = hp.nside2pixarea(nside)
    pix_deg2 = hp.nside2pixarea(nside, degrees=True)
    area_sr = n_covered * pix_sr
    area_deg2 = n_covered * pix_deg2
    sky_fraction = n_covered / npix if npix else 0.0

    return {
        "mission": mission.upper(),
        "path": map_path,
        "nside": nside,
        "npix": npix,
        "n_covered": n_covered,
        "area_deg2": area_deg2,
        "area_sr": area_sr,
        "sky_fraction": sky_fraction,
        "total_exptime": float(np.sum(hp_map)),
        "min_exptime": min_exptime,
    }


def format_result(result: dict, show_pixels: bool = False) -> str:
    """Format coverage statistics for a single mission."""
    lines = [
        (
            f"{result['mission']} coverage "
            f"({result['path']}, nside={result['nside']})"
        ),
    ]
    if show_pixels:
        lines.append(
            f"  covered pixels : {result['n_covered']:,} / {result['npix']:,}"
        )
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
    return "\n".join(lines)


def format_all_results(results: list[dict], show_pixels: bool = False) -> str:
    """Format coverage statistics as a table for multiple missions."""
    if show_pixels:
        header = (
            f"{'Mission':<8} {'Pixels':>12} {'deg^2':>12} "
            f"{'sr':>10} {'Fraction':>10} {'Percent':>9}"
        )
    else:
        header = (
            f"{'Mission':<8} {'deg^2':>12} "
            f"{'sr':>10} {'Fraction':>10} {'Percent':>9}"
        )
    separator = "-" * len(header)
    rows = [header, separator]

    for result in results:
        if show_pixels:
            rows.append(
                f"{result['mission']:<8} "
                f"{result['n_covered']:>12,} "
                f"{result['area_deg2']:>12.2f} "
                f"{result['area_sr']:>10.4f} "
                f"{result['sky_fraction']:>10.5f} "
                f"{result['sky_fraction'] * 100:>8.3f}%"
            )
        else:
            rows.append(
                f"{result['mission']:<8} "
                f"{result['area_deg2']:>12.2f} "
                f"{result['area_sr']:>10.4f} "
                f"{result['sky_fraction']:>10.5f} "
                f"{result['sky_fraction'] * 100:>8.3f}%"
            )

    return "\n".join(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Report sky coverage for a MAST mission from cached HEALPix maps."
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
        help="Report coverage for every mission with a cached map",
    )
    parser.add_argument(
        "--data-dir",
        default="data",
        help="Directory containing <mission>_map.fits files (default: data)",
    )
    parser.add_argument(
        "--min-exptime",
        type=float,
        default=0.0,
        help="Minimum exposure time for a pixel to count as covered",
    )
    parser.add_argument(
        "--pixels",
        action="store_true",
        help="Include covered pixel counts in the output",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.all and args.mission:
        parser.error("Specify either a mission name or --all, not both.")
    if not args.all and not args.mission:
        parser.error("Specify a mission name or use --all.")

    if args.all:
        missions = available_missions(args.data_dir)
        if not missions:
            parser.error(f"No mission maps found in {args.data_dir}")
        results = [
            mission_area(mission, data_dir=args.data_dir, min_exptime=args.min_exptime)
            for mission in missions
        ]
        print(format_all_results(results, show_pixels=args.pixels))
        return 0

    result = mission_area(
        args.mission,
        data_dir=args.data_dir,
        min_exptime=args.min_exptime,
    )
    print(format_result(result, show_pixels=args.pixels))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

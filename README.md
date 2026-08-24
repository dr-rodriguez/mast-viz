# mast-viz

Coverage map visualizations of MAST data.

Uses healpy and cdshealpix to generate exposure maps and astronify to sonify them.

# Setup

Use uv to install the dependencies.

```bash
uv sync
```

Copy the .env.example file to .env and update the environment variables.

```bash
cp .env.example .env
```

# Updating Data

Run individual lines in the fetch_data.py file to update the data.   
Pending work is to split this into scripts for each mission and organize/document it better.

# Coverage Area

Report sky coverage for a mission from its cached HEALPix map in `data/<mission>_map.fits`.
Coverage is reported in square degrees, steradians, and fraction of the full sky.

```bash
uv run python -m mast_viz.mission_area hst
uv run python -m mast_viz.mission_area --all
uv run python -m mast_viz.mission_area jwst --min-exptime 1000
```

For more accurate coverage, use `mission_area_exact`, which unions the exact
`s_region` footprints from `data/<mission>.h5` into a MOC. This is slower on
the first run but caches the result to `data/<mission>_moc_o<order>.fits`.

Default HEALPix orders per mission: hst/jwst 18, sdss/galex/hlsp 16, ps1/kepler
14, tess 11. Override with `--order`. Use `--compare` to print the map-based
estimate alongside the exact value.

```bash
uv run python -m mast_viz.mission_area_exact hst
uv run python -m mast_viz.mission_area_exact --all
uv run python -m mast_viz.mission_area_exact hst --compare
uv run python -m mast_viz.mission_area_exact tess --order 13 --recompute
```

The map-based script can overestimate coverage (e.g. HST map-based ~2216 deg^2
vs exact ~77 deg^2) because nside=256 pixels that merely intersect a footprint
are counted at ~0.0525 deg^2 each.

# Generate Frames

Use scripts like galex_movie.py to generate the frames for each mission.

# Generate Movies

Refer to the script movie/movie_script.sh to see how to call ffmpeg to generate the movies.

# Generate Sonifications

Refer to mast_sound.py to see how to generate the sonifications.
This will require the astronify package and extra setup. See https://astronify.readthedocs.io/en/latest/astronify/install.html

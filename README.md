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

Run individual lines in the script.py file to update the data.   
Pending work is to split this into scripts for each mission and organize/document it better.

# Generate Frames

Use scripts like galex_movie.py to generate the frames for each mission.

# Generate Movies

Refer to the script movie/movie_script.sh to see how to call ffmpeg to generate the movies.

# Generate Sonifications

Refer to mast_sound.py to see how to generate the sonifications.
# Use astronify to make sound of the MAST data
import matplotlib
matplotlib.use('Qt5Agg')  # avoids crashing MacOS Mojave
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy.time import Time
from astropy.table import Table
from astronify.series import SoniSeries
plt.interactive(False)

# Read the data
df = pd.read_csv('data/mast_time.csv')
t = Table.from_pandas(df)
t['year'] = [x.datetime.strftime('%Y-%m') for x in Time(t['date'], format='mjd')]

# Quick plots
plt.plot(t['week'], t['area'])  # area on sky over time
plt.plot(t['week'], t['obs_counts'])  # number of observations over time
plt.plot(t['week'], t['exp_counts'])  # exposure time over time

plt.plot(t['year'], t['area'])  # area on sky over time
plt.xticks(rotation=90)
ax = plt.gca()
ax.set_xticks(ax.get_xticks()[::10])
plt.tight_layout()
plt.savefig('image/mast_area.png', dpi=300)
plt.close()

plt.plot(t['year'], np.log10(t['area']))  # area on sky over time
plt.xticks(rotation=90)
ax = plt.gca()
ax.set_xticks(ax.get_xticks()[::10])
plt.tight_layout()
plt.savefig('image/mast_logarea.png', dpi=300)
plt.close()

plt.plot(t['year'], t['exp_counts'])  # exposure time over time
plt.xticks(rotation=90)
ax = plt.gca()
ax.set_xticks(ax.get_xticks()[::10])
plt.tight_layout()
plt.savefig('image/mast_exptime.png', dpi=300)
plt.close()

plt.plot(t['year'], t['obs_counts'])  # observation counts over time
plt.xticks(rotation=90)
ax = plt.gca()
ax.set_xticks(ax.get_xticks()[::10])
plt.tight_layout()
plt.savefig('image/mast_counts.png', dpi=300)
plt.close()


# Notable dates: 1997 for STIS and NICMOS, 2002 for ACS, 2009 for WFC3

# Sound
frame_rate = 20.  # note spacing will be 1/frame_rate
sonidf = SoniSeries(t, time_col='week', val_col='area')
print(sonidf.pitch_mapper.pitch_map_args)  # check parameters
sonidf.note_spacing = 1/frame_rate  # median seconds between notes; default is 0.01; my video uses 20 frames/sec
sonidf.pitch_mapper.pitch_map_args["stretch"] = "log"  # or linear
# sonidf.pitch_mapper.pitch_map_args["center_pitch"] = 800  # overall higher pitch
sonidf.sonify()
sonidf.play()
sonidf.stop()
sonidf.write("movie/mast_logarea.wav")

sonidf = SoniSeries(t, time_col='week', val_col='obs_counts')
sonidf.note_spacing = 1/frame_rate
sonidf.sonify()
sonidf.write("movie/mast_counts.wav")

sonidf = SoniSeries(t, time_col='week', val_col='exp_counts')
# sonidf.pitch_mapper.pitch_map_args["stretch"] = "log"
sonidf.note_spacing = 1/frame_rate
sonidf.sonify()
sonidf.write("movie/mast_exptime.wav")


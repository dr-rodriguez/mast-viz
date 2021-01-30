# Scripts to generate quick movies with ffmpeg

# GALEX
~/Software/ffmpeg -r 10 -i ./galex/galex_frame%06d.png -vcodec mpeg4 -q:v 3 -y galex.mp4

# MAST
~/Software/ffmpeg -r 20 -i ./mast/mast_frame%06d.png -vcodec mpeg4 -q:v 3 -y mast.mp4

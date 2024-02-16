### CONFIG FILE ###

import os


#~ PATH DEFINITIONS
DIRNAME = os.path.dirname(__file__)
DATAPATH = os.path.join("/home/guillaumechagnaud/TRAVAIL/POSTDOC/IGE-Secheresse","data")
SHAPEDIR = os.path.join(DATAPATH,"shapes")
OUTDIR = os.path.join(DIRNAME,"output")
DATADIR = os.path.join(OUTDIR,"data")
FIGDIR = os.path.join(OUTDIR,"figures")



# SAFRAN PROPERTIES

res = 8000 # horizontal resolution (m)


#~ FRANCE TERRITORY BOUNDARIES (m)

x_min = 60000
x_max = 1196000
x_size = x_max - x_min
y_min = 1617000
y_max = 2681000
y_size = y_max - y_min


# Restricted boundaries -> sub-domains contain enough data
lat_min = 43.
lat_max = 50.
lon_min = -2.
lon_max = 7.5

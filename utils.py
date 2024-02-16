"""Some functions"""

from pyproj import Transformer
from collections import Iterable


### CST ###

# SAFRAN MAP PROJECTION
# grid_mapping_name:              lambert_conformal_conic
# standard_parallel:              [45.89892 47.69601]
# longitude_of_central_meridian:  2.337229
# latitude_of_projection_origin:  46.8
# false_easting:                  600000.0
# false_northing:                 2200000.0
proj = "+proj=lcc +lat_0=46.8 +lon_0=2.337229 +lat_1=45.89892 +lat_2=47.69601 +x_0=600000. +y_0=2200000. +units=m +datum=NAD83 +no_defs"


### FUNCTIONS ###

def transform_latlon_to_lcc(lat, lon):
    """Transform lat, lon coordinates into Lambert Conic Conformal with the SAFRAN proj"""
    transformer_latlon_to_lcc = Transformer.from_proj("EPSG:4326", proj)
    out = transformer_latlon_to_lcc.transform(lat, lon)  # return lon, lat in LCC coordinates
    out = (out[1], out[0])  # return lat, lon

    return(out)


def transform_lcc_to_latlon(x, y):
    """Transform Lambert Conic Conformal with the SAFRAN proj coordinates into lat, lon"""
    transformer_lcc_to_latlon = Transformer.from_proj(proj, "EPSG:4326")
    out = transformer_lcc_to_latlon.transform(x, y)  # return lon, lat in WGS84
    out = (out[1], out[0])  # return lat, lon

    return(out)


def transform_lcc_to_lambert93(x, y):
    """Transform Lambert Conic Conformal with the SAFRAN proj coordinates into Lambert 2"""
    transformer_lcc_to_lamb2 = Transformer.from_proj(proj, "EPSG:2154")
    out = transformer_lcc_to_lamb2.transform(x, y)  # return x, y in Lambert 2
    out = (out[0], out[1])

    return(out)


def rgb_to_hex(rgb):
    """Transform color rgb values to hex"""
    if isinstance(rgb, dict):
        r = rgb['r']
        g = rgb['g']
        b = rgb['b']
    
    elif isinstance(rgb, Iterable):
        if len(rgb) == 3:
            pass
        elif len(rgb) == 4:
            assert rgb[3] == 256 - 1, "alpha paramter is expected to be 1, got {0}".format(rgb)
        else:
            raise ValueError("3 colors are expected, got {0}".format(rgb))

        r = rgb[0]
        g = rgb[1]
        b = rgb[2]
            
    else:
        raise ValueError("Not implementd")
        
    r = int(round(r))
    b = int(round(b))
    g = int(round(g))
    res = "#{0:02x}{1:02x}{2:02x}".format(r, g, b)
    assert len(res) == 7, "expected hexa, got {0}".format(res)

    return res



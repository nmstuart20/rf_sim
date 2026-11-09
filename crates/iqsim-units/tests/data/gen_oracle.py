"""Generate MGRS / UTM reference data for iqsim-units tests.

Usage: python -I gen_oracle.py <out_dir>

1. mgrs_geoconvert.txt: MGRS strings (encoded by the NGA GEOTRANS-based `mgrs`
   package from chosen lat/lon points) decoded by GeographicLib's published
   GeoConvert web service (square center lat/lon and UTM easting/northing).
2. utm_pyproj.txt: random UTM points inverted to lat/lon with PROJ (pyproj).
"""

import html
import random
import re
import sys
import time
import urllib.parse
import urllib.request

import mgrs
import pyproj

out = sys.argv[1]
rng = random.Random(20261108)
m = mgrs.MGRS()

# Hand-picked places: exception zones, both hemispheres, band edges, antimeridian.
fixed = [
    (38.8895, -77.0352, 5),    # Washington DC
    (21.3069, -157.8583, 5),   # Honolulu
    (-33.8568, 151.2153, 5),   # Sydney
    (-54.8019, -68.3030, 4),   # Ushuaia, band F
    (-79.5, 10.0, 3),          # band C
    (83.6, -32.66, 5),         # band X, Greenland
    (60.5, 5.5, 5),            # zone 32V (widened)
    (61.0, 2.5, 5),            # 31V (narrowed)
    (75.0, 8.0, 5),            # 31X (widened)
    (75.0, 10.0, 5),           # 33X
    (78.0, 22.0, 5),           # 35X
    (79.0, 40.0, 5),           # 37X
    (0.0001, 0.0001, 5),       # just north of the equator
    (-0.0001, 0.0001, 5),      # just south of the equator
    (0.0, 3.0, 5),             # zone edge at the equator
    (40.0, 179.9, 5),          # antimeridian, zone 60
    (40.0, -179.9, 5),         # zone 1
    (51.4779, -0.0015, 1),     # Greenwich, 10 km
    (51.4779, -0.0015, 0),     # Greenwich, 100 km square
    (-26.2041, 28.0473, 2),    # Johannesburg
    (64.1466, -21.9426, 5),    # Reykjavik
    (-77.85, 166.67, 5),       # McMurdo, band C
    (8.0, 81.0, 5),            # band N bottom edge region
    (71.9, -100.0, 5),         # band W top
    (72.1, -100.0, 5),         # band X bottom
    (56.1, 3.1, 5),            # 32V lower-left corner region
]
for _ in range(34):
    lat = rng.uniform(-80.0, 84.0)
    lon = rng.uniform(-180.0, 180.0)
    fixed.append((lat, lon, rng.choice([0, 1, 2, 3, 4, 5, 5, 5])))

PAT = re.compile(
    r"lat lon = .*?>(-?[0-9.]+) (-?[0-9.]+)</a>.*?"
    r"UPS</a> = <font color=\"Blue\">(\d+)([ns]) ([0-9.]+) ([0-9.]+)</font>",
    re.S,
)

rows = []
for lat, lon, prec in fixed:
    ref = m.toMGRS(lat, lon, MGRSPrecision=prec)
    if isinstance(ref, bytes):
        ref = ref.decode()
    q = urllib.parse.urlencode(
        {"input": ref, "option": "Submit", "format": "g", "prec": "6", "zone": "-3"}
    )
    url = "https://geographiclib.sourceforge.io/cgi-bin/GeoConvert?" + q
    body = urllib.request.urlopen(url, timeout=30).read().decode()
    hit = PAT.search(body)
    if not hit:
        raise SystemExit(f"no result for {ref}")
    clat, clon, zone, hemi, e, n = hit.groups()
    rows.append(f"{ref} {clat} {clon} {zone}{hemi} {e} {n}")
    time.sleep(0.3)

with open(f"{out}/mgrs_geoconvert.txt", "w") as f:
    f.write(
        "# MGRS reference points decoded by GeographicLib GeoConvert 2.5.2\n"
        "# (https://geographiclib.sourceforge.io/cgi-bin/GeoConvert, fetched 2026-11-08).\n"
        "# MGRS strings were encoded from chosen lat/lon points with the Python `mgrs`\n"
        "# package (NGA GEOTRANS). GeoConvert returns the CENTER of the grid square.\n"
        "# Columns: mgrs center_lat_deg center_lon_deg zone+hemisphere center_easting_m center_northing_m\n"
    )
    f.write("\n".join(rows) + "\n")

# UTM inverse against PROJ.
lines = []
for _ in range(400):
    zone = rng.randint(1, 60)
    south = rng.random() < 0.5
    e = rng.uniform(150_000.0, 850_000.0)
    n = rng.uniform(1_200_000.0, 9_999_000.0) if south else rng.uniform(0.0, 9_300_000.0)
    # Round first so the file holds exactly the inputs PROJ inverted.
    e = float(f"{e:.4f}")
    n = float(f"{n:.4f}")
    crs = pyproj.CRS.from_dict({"proj": "utm", "zone": zone, "south": south, "ellps": "WGS84"})
    t = pyproj.Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    lon, lat = t.transform(e, n)
    lines.append(f"{zone}{'s' if south else 'n'} {e:.4f} {n:.4f} {lat:.12f} {lon:.12f}")
with open(f"{out}/utm_pyproj.txt", "w") as f:
    f.write(
        f"# UTM (WGS84) inverse reference points from PROJ {pyproj.proj_version_str} via pyproj {pyproj.__version__}.\n"
        "# Columns: zone+hemisphere easting_m northing_m lat_deg lon_deg\n"
    )
    f.write("\n".join(lines) + "\n")
print(len(rows), len(lines))

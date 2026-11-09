"""Generate WGS84 geodetic <-> ECEF reference data for iqsim-geo tests.

Usage: python -I gen_oracle.py <out_dir>

Two independent implementations compute every row:
- PROJ (through pyproj), EPSG:4979 (WGS84 3D geographic) <-> EPSG:4978 (geocentric).
- An mpmath reference at 60 significant digits: the closed-form forward map, and
  for the inverse the fixed-point iteration phi <- atan2(z + e^2 N sin(phi), p),
  run 80 times (contraction factor <= e^2 ~ 0.0067 outside the evolute, so the
  result is exact to far below f64 resolution). Inputs are the exact values of the
  f64s in the file; outputs are rounded once to f64.

Files (one row per line, whitespace-separated, shortest round-trip floats):
1. geodetic_forward.txt: lat_deg lon_deg h_m | x y z (PROJ) | x y z (mpmath)
2. geodetic_inverse.txt: x y z | lat lon h (PROJ) | lat lon h (mpmath)

The points are hand-picked edge cases (poles, equator, the antimeridian on both
sides, ocean-trench depth, airliner and LEO/GEO heights, a lunar distance) plus
seeded random points. Not run in CI; the outputs are checked in.
"""

import math
import random
import sys

import mpmath
import pyproj
from mpmath import mp, mpf

out = sys.argv[1]
rng = random.Random(20261109)

mp.dps = 60
A = mpf(6378137)
F = 1 / (mpf("298.257223563"))
E2 = F * (2 - F)
DEG = mp.pi / 180


def mp_forward(lat, lon, h):
    phi, lam, h = mpf(lat) * DEG, mpf(lon) * DEG, mpf(h)
    n = A / mp.sqrt(1 - E2 * mp.sin(phi) ** 2)
    return ((n + h) * mp.cos(phi) * mp.cos(lam),
            (n + h) * mp.cos(phi) * mp.sin(lam),
            (n * (1 - E2) + h) * mp.sin(phi))


def mp_inverse(x, y, z):
    x, y, z = mpf(x), mpf(y), mpf(z)
    p = mp.sqrt(x * x + y * y)
    phi = mp.atan2(z, p)
    for _ in range(80):
        n = A / mp.sqrt(1 - E2 * mp.sin(phi) ** 2)
        phi = mp.atan2(z + E2 * n * mp.sin(phi), p)
    h = p * mp.cos(phi) + z * mp.sin(phi) - A * mp.sqrt(1 - E2 * mp.sin(phi) ** 2)
    lam = mp.atan2(y, x) if p != 0 else mpf(0)
    return phi / DEG, lam / DEG, h


fwd = pyproj.Transformer.from_crs("EPSG:4979", "EPSG:4978", always_xy=True)
inv = pyproj.Transformer.from_crs("EPSG:4978", "EPSG:4979", always_xy=True)

heights = [-11_000.0, -100.0, 0.0, 10.0, 10_668.0, 400_000.0, 2_000_000.0,
           35_786_000.0, 384_400_000.0]
fixed = []
for h in heights:
    for lat, lon in [
        (0.0, 0.0), (0.0, 90.0), (0.0, -90.0), (0.0, 180.0), (0.0, -180.0),
        (90.0, 0.0), (-90.0, 0.0), (89.999999, 45.0), (-89.999999, -135.0),
        (45.0, 179.9999999), (45.0, -179.9999999), (38.8895, -77.0352),
        (-33.8568, 151.2153), (71.0, 25.0), (1e-9, 1e-9), (-1e-9, -1e-9),
    ]:
        fixed.append((lat, lon, h))

rand = []
for _ in range(400):
    lat = math.degrees(math.asin(rng.uniform(-1.0, 1.0)))
    lon = rng.uniform(-180.0, 180.0)
    band = rng.random()
    if band < 0.5:
        h = rng.uniform(-11_000.0, 20_000.0)
    elif band < 0.85:
        h = rng.uniform(20_000.0, 2_000_000.0)
    else:
        h = rng.uniform(2_000_000.0, 400_000_000.0)
    rand.append((lat, lon, h))


def r(v):
    return repr(float(v))


header = (
    f"# pyproj {pyproj.__version__} (PROJ {pyproj.proj_version_str}), "
    f"mpmath {mpmath.__version__} at {mp.dps} digits, seed 20261109\n"
)

with open(f"{out}/geodetic_forward.txt", "w") as f:
    f.write(header)
    f.write("# lat_deg lon_deg h_m  x y z (PROJ)  x y z (mpmath)\n")
    for lat, lon, h in fixed + rand:
        px, py, pz = fwd.transform(lon, lat, h)
        gx, gy, gz = (float(v) for v in mp_forward(lat, lon, h))
        f.write(" ".join(r(v) for v in (lat, lon, h, px, py, pz, gx, gy, gz)) + "\n")

with open(f"{out}/geodetic_inverse.txt", "w") as f:
    f.write(header)
    f.write("# x y z  lat lon h (PROJ)  lat lon h (mpmath)\n")
    pts = []
    for lat, lon, h in fixed + rand:
        pts.append(tuple(float(v) for v in mp_forward(lat, lon, h)))
    # Random ECEF points by direction and radius, including near the axis.
    for _ in range(200):
        x, y, z = (rng.gauss(0, 1) for _ in range(3))
        n = math.sqrt(x * x + y * y + z * z)
        rad = rng.choice([rng.uniform(6.34e6, 6.40e6), rng.uniform(6.4e6, 4.3e7)])
        pts.append((x / n * rad, y / n * rad, z / n * rad))
    for _ in range(20):
        pts.append((rng.uniform(-1e-3, 1e-3), rng.uniform(-1e-3, 1e-3),
                    rng.choice([-1, 1]) * rng.uniform(6.3e6, 4.3e7)))
    for x, y, z in pts:
        plon, plat, ph = inv.transform(x, y, z)
        gl = [float(v) for v in mp_inverse(x, y, z)]
        f.write(" ".join(r(v) for v in (x, y, z, plat, plon, ph, gl[0], gl[1], gl[2])) + "\n")

# 3. eci_angles.txt: t_ps frac theta_rad -- the Earth rotation angle
#    omega_E * (t_ps + frac*2^-64) * 1e-12 reduced to [0, 2*pi), omega_E =
#    7292115e-11 rad/s exactly, computed in mpmath and rounded once to f64.
with open(f"{out}/eci_angles.txt", "w") as f:
    f.write(header)
    f.write("# t_ps frac theta_rad\n")
    omega = mpf(7292115) * mpf(10) ** -11
    two_pi = 2 * mp.pi
    ts = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, 2**64 - 1), (86_164_090_500_000, 0),
          (86_400 * 10**12, 0), (100 * 86_400 * 10**12, 12345678901234567890),
          (-100 * 86_400 * 10**12, 0), (2**63 - 1, 2**64 - 1), (-(2**63), 0),
          (5858921581365495166, 0)]  # review: 1.52 ulp with the old conversion
    for _ in range(200):
        ts.append((rng.randrange(-(2**63), 2**63), rng.randrange(0, 2**64)))
    for t, fr in ts:
        secs = (mpf(t) + mpf(fr) / mpf(2) ** 64) * mpf(10) ** -12
        th = mp.fmod(omega * secs, two_pi)
        if th < 0:
            th += two_pi
        f.write(f"{t} {fr} {r(th)}\n")

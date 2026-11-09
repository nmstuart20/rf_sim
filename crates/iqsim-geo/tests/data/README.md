# iqsim-geo reference data

| File | Source | Used by |
| --- | --- | --- |
| `airsport_modecascii.txt` | AirSport Avionics, Darryl Phillips, "Table of Mode A and Mode C encoder data" (1998, archived 26 July 2012): <https://web.archive.org/web/20120726003224/http://www.airsport-corp.com/modecascii.txt>, cited by Wikipedia's "Gillham code" article. Unmodified. SHA-256 `b750d7d1a717c4e9faf1fa8897450ec440e219101c634288326827cb4de88c21`. 1280 rows, −1200 ft to 126 700 ft. | `tests/gillham.rs` |
| `geodetic_forward.txt`, `geodetic_inverse.txt` | `gen_oracle.py`: PROJ 9.8.1 (pyproj 3.8.0, EPSG:4979 ↔ EPSG:4978) and a 60-digit mpmath 1.4.1 reference. | `tests/geodesy.rs` |
| `eci_angles.txt` | `gen_oracle.py`: Earth rotation angle in mpmath. | `tests/geodesy.rs` |

Regenerate the generated files with `python -I gen_oracle.py .` in a virtual environment with `pyproj` and `mpmath`. The script is seeded, so the output is reproducible for the same tool versions. It is not run in CI.

PROJ's geocentric inverse is a fixed-step approximation (about 10⁻⁵ m at the surface, centimetres at 2000 km), so the tests hold the inverse to the mpmath reference and use PROJ only as a loose cross-check near the surface. The forward map agrees between PROJ and mpmath to 8·10⁻⁸ m.

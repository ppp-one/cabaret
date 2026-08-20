# Local catalog

`GaiaQuery` supports querying offline SQLite catalogs in addition to online
TAP services. While TAP services offer access to the full Gaia archive, local SQLite catalogs provide significantly lower latency and faster query execution,
making them ideal for large-scale simulations or working without an internet connection.

## Using a local SQLite catalog

Create a `GaiaSQLiteSource` to point Cabaret at a local SQLite file and
optional table name. The table is expected to expose Gaia-like columns named
similarly to the TAP service (for example: `ra`, `dec`, `phot_g_mean_mag`,
`h_m`, ...). Example:

```python
import cabaret
from astropy.coordinates import SkyCoord

sqlite_source = cabaret.GaiaSQLiteSource(
    database="/data/catalogs/gaia_subset.sqlite",
    table="gaia_sources",
)

center = SkyCoord(ra=323.36, dec=-0.82, unit="deg")
table = cabaret.GaiaQuery.query(
    center=center,
    radius=0.1,  # degrees
    filter_bands=[cabaret.Filters.G, cabaret.Filters.J],
    tap_source=sqlite_source,  # or "sqlite:///path/to/file.sqlite"
)
```

You may pass the `tap_source` as either a `GaiaSQLiteSource` instance or a
string. Accepted string forms include a full SQLite URI (`sqlite:///...`) or a
plain path ending in `.db`, `.sqlite` or `.sqlite3`.

### Changing the default TAP source

You can change the module-wide default TAP source so that calls which omit
`tap_source` use your local SQLite file by default. For example:

```python
import cabaret

cabaret.GaiaQuery.DEFAULT_TAP_SOURCE = cabaret.GaiaSQLiteSource(
    database="/data/catalogs/gaia_subset.sqlite",
    table="table_name",
)

# Subsequent calls will use the SQLite database by default
image = cabaret.Observatory().generate_image(
    ra=12.33230,  # right ascension in degrees
    dec=30.4343,  # declination in degrees
    exp_time=10,  # exposure time in seconds
)
```

## Sharded (declination-ring) catalogs

Some offline catalogs are split into declination "rings" (shards) named like
`-90_-89`, `-89_-88`, etc. If the configured `table` name is not present in
the database, Cabaret will auto-detect these ring tables and only query the
shards that overlap your requested declination range. This keeps queries fast
when only a small declination slice is required.

Prebuilt sharded SQLite databases of the Gaia-2MASS crossmatch catalogs containing G and J magnitudes used at the SPECULOOS observatory can be downloaded here: [Gaia-2MASS Local Sharded SQLite Catalogue](https://zenodo.org/records/18214672). In the example below, we download the `gaia_tmass_7_jm_cut.db` catalog (~11 MB) and query a region across `RA=0`.

```python
import requests

import cabaret

url = "https://zenodo.org/records/18214672/files/gaia_tmass_7_jm_cut.db"
local_path = "gaia_tmass_7_jm_cut.db"

with open(local_path, "wb") as f:
    f.write(requests.get(url).content)

table = cabaret.GaiaQuery.query(
    bounds=(359.5, 0.5, -1.0, 1.0),  # ra_min, ra_max wrap across 0 deg
    filter_bands=["G", "J"],
    tap_source=local_path,
)
```

## RA wraparound

RA values are interpreted modulo 360. If `ra_min <= ra_max` the interval is
treated as contiguous; if `ra_min > ra_max` the interval is interpreted as
wrapping across RA=0 (for example `(350.0, 10.0, -5.0, 5.0)`). Both the
TAP and SQLite query code handle wraparound automatically, selecting rows
whose RA falls into the wrapped interval.

## See also

- {doc}`GaiaQuery API Reference <generated/cabaret.GaiaQuery>`
- Prebuilt sharded Gaia+2MASS SQLite catalogs: [ppp-one/gaia-tmass-sqlite](https://github.com/ppp-one/gaia-tmass-sqlite)

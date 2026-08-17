# exoplanet-hunt

Searching for new transiting exoplanet candidates in the freshest public TESS
data — stars that are in **no TOI list and no confirmed-planet table**.

## The idea

TESS keeps observing, but the official TOI alert pipeline lags the newest
sector releases, and small faint stars get the least follow-up scrutiny.
This project searches the most recently released TESS sector (currently
**sector 104**, observed May 2026, 12,996 two-minute-cadence light curves)
and cross-matches every signal against the live ExoFOP TOI catalog and the
NASA Exoplanet Archive confirmed-planet table. A strong transit signal on a
star in neither list is a genuinely new candidate.

Target selection favors M/K dwarfs (R* < 0.75 R_sun, Teff < 4800 K,
Tmag < 14): the smaller the star, the deeper the transit a given planet
produces.

## Pipeline

1. `build_sample.py` — parse the sector's MAST bulk-download manifest, drop
   every TIC already flagged as a TOI or confirmed host, pull stellar
   parameters from the TESS Input Catalog, rank the remainder
   (12,689 unflagged stars → 400-target M/K-dwarf sample).
2. `hunt.py` — per star: download the SPOC PDCSAP light curve, clean,
   flatten, Box Least Squares period search (0.5–13 d), then vet:
   - detection strength (SDE ≥ 8) and ≥ 2 distinct transits covered
   - depth ≤ 10 % (deeper ⇒ eclipsing binary)
   - odd/even transit depth consistency (binaries alternate; planets don't)
   - no significant secondary eclipse at phase 0.5
   Survivors land in `candidates/candidates.csv` with a 4-panel diagnostic
   figure in `figures/`. The full sweep record is `data/results.csv`.
3. `validation/find_exoplanets.py` — the same method run on known planets.
   It recovers WASP-18 b (P = 0.9417 d vs 0.9415 d, S/N 25) and
   LHS 3844 b (P = 0.4629 d, exact, S/N 35). It does **not** recover the
   ~300 ppm transit of π Men c from one sector — the pipeline is sensitive
   to transits ≳ 1000 ppm, i.e. roughly Earth-size and larger planets
   around the M dwarfs targeted here.

## Reproduce

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
# fetch catalogs (not committed; live versions change daily):
curl -s "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv" -o data/toi.csv
curl -s "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select+pl_name,hostname,tic_id,disc_facility,disc_year+from+pscomppars&format=csv" -o data/confirmed.csv
curl -s "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_104_lc.sh" -o data/sector104_lc.sh
.venv/bin/python build_sample.py
.venv/bin/python hunt.py     # resumable; skips already-searched stars
```

## Honest caveats

- A surviving signal is a **candidate**, not a discovery. Ruling out
  background eclipsing binaries needs pixel-level centroid checks and
  usually ground-based follow-up (that's what ExoFOP is for).
- The SPOC pipeline team searches these same light curves; the edge here is
  timing (TOI alerts lag new sectors) and depth of scrutiny on faint cool
  dwarfs, not access to secret data.
- Single-sector search: only periods ≲ 13 d give the required two transits.

Data: NASA TESS mission via MAST; TOI catalog via ExoFOP-TESS; confirmed
planets via the NASA Exoplanet Archive.

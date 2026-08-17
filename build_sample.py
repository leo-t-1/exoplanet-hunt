"""Build the under-searched target sample from TESS sector 104.

Parses TIC IDs from the MAST bulk-download manifest, removes every star that
is already a TOI or a confirmed planet host, pulls stellar parameters from
the TESS Input Catalog, and ranks the remainder to favor small cool stars
(deepest transits for a given planet size).

Output: data/targets.csv  (rank-ordered, with download URLs)
Run:    .venv/bin/python build_sample.py
"""

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astroquery.mast import Catalogs

warnings.filterwarnings("ignore")

DATA = Path(__file__).parent / "data"
N_KEEP = 400  # how many ranked targets to write out

# --- parse sector manifest ------------------------------------------------
manifest = (DATA / "sector104_lc.sh").read_text()
pat = re.compile(r"-o (tess\d+-s0104-(\d{16})-\d+-s_lc\.fits) (https://\S+)")
rows = [(int(tic), fname, url) for fname, tic, url in pat.findall(manifest)]
sector = pd.DataFrame(rows, columns=["tic", "fname", "url"]).drop_duplicates("tic")
print(f"sector 104 light curves : {len(sector)}")

# --- exclusion sets: anything humanity has already flagged ----------------
toi = pd.read_csv(DATA / "toi.csv", usecols=["TIC ID"])
toi_tics = set(toi["TIC ID"].astype(int))
conf = pd.read_csv(DATA / "confirmed.csv")
conf_tics = set(
    int(m.group(1))
    for s in conf["tic_id"].dropna()
    if (m := re.search(r"(\d+)", str(s)))
)
known = toi_tics | conf_tics
print(f"known TOIs              : {len(toi_tics)}")
print(f"confirmed-host TICs     : {len(conf_tics)}")

fresh = sector[~sector["tic"].isin(known)].reset_index(drop=True)
print(f"unflagged sector targets: {len(fresh)}")

# --- stellar parameters from the TESS Input Catalog (chunked) -------------
params = []
tics = fresh["tic"].tolist()
for i in range(0, len(tics), 1000):
    chunk = [str(t) for t in tics[i : i + 1000]]
    tab = Catalogs.query_criteria(catalog="TIC", ID=chunk)
    params.append(tab.to_pandas()[["ID", "Tmag", "Teff", "rad", "mass", "d"]])
    print(f"  TIC query {i + len(chunk)}/{len(tics)}")
tic_params = pd.concat(params)
tic_params["ID"] = tic_params["ID"].astype(int)
df = fresh.merge(tic_params, left_on="tic", right_on="ID", how="left")

# --- rank: small cool dwarfs first, then bright ---------------------------
ok = df[(df["rad"] > 0.08) & (df["rad"] < 0.75) & (df["Teff"] < 4800) & (df["Tmag"] < 14)]
ok = ok.sort_values(["rad", "Tmag"]).head(N_KEEP)
print(f"ranked M/K-dwarf sample : {len(ok)}")

out = DATA / "targets.csv"
ok[["tic", "Tmag", "Teff", "rad", "mass", "d", "fname", "url"]].to_csv(out, index=False)
print(f"wrote {out}")

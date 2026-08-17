"""Second-stage vetting of BLS candidates from hunt.py.

Physical + statistical cuts on data/results.csv, then TLS refinement of the
survivors (limb-darkened transit model, independent SDE + false-alarm
probability, odd/even test).

Cuts:
  - duty cycle: duration/period must be < 0.10 (a real transit around a
    small star occupies a few percent of the orbit; long "transits" at short
    periods are detrending residuals / scattered light)
  - period off the search-grid edge (> 0.52 d)
  - no strong phase-0.5 anomaly of either sign (|secondary_sigma| < 5)
  - systematics dedup: if several stars share a period within 2% and an
    epoch within 0.1 d, that's the spacecraft, not planets — drop them all

Output: candidates/vetted.csv
Run:    .venv/bin/python vet.py
NB: body lives under a __main__ guard — TLS spawns multiprocessing workers
that re-import this module on macOS.
"""

import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import lightkurve as lk
from transitleastsquares import transitleastsquares

warnings.filterwarnings("ignore")

ROOT = Path(__file__).parent
DATA, CAND = ROOT / "data", ROOT / "candidates"
CACHE = DATA / "lc_cache"


def main():
    CACHE.mkdir(exist_ok=True)
    res = pd.read_csv(DATA / "results.csv")
    cand = res[res["verdict"] == "CANDIDATE"].copy()
    print(f"raw BLS candidates      : {len(cand)}")

    cand["duty"] = cand["duration_hr"] / 24.0 / cand["period"]
    cand = cand[cand["duty"] < 0.10]
    print(f"after duty-cycle cut    : {len(cand)}")
    cand = cand[cand["period"] > 0.52]
    print(f"after grid-edge cut     : {len(cand)}")
    cand = cand[cand["secondary_sigma"].abs() < 5]
    print(f"after phase-0.5 cut     : {len(cand)}")

    drop = set()
    c = cand.reset_index(drop=True)
    for i in range(len(c)):
        for j in range(i + 1, len(c)):
            pi, pj = c.loc[i, "period"], c.loc[j, "period"]
            if abs(pi - pj) / pi < 0.02 and abs(c.loc[i, "t0"] - c.loc[j, "t0"]) < 0.1:
                drop |= {c.loc[i, "tic"], c.loc[j, "tic"]}
    cand = cand[~cand["tic"].isin(drop)]
    print(f"after systematics dedup : {len(cand)} (dropped {len(drop)} shared-signal stars)")

    targets = pd.read_csv(DATA / "targets.csv").set_index("tic")
    rows = []
    for r in cand.itertuples():
        tic = int(r.tic)
        tinfo = targets.loc[tic]
        fpath = CACHE / tinfo["fname"]
        if not fpath.exists():
            urllib.request.urlretrieve(tinfo["url"], fpath)
        lc = lk.read(fpath, flux_column="pdcsap_flux")
        lc = lc.remove_nans().normalize().remove_outliers(sigma_upper=3, sigma_lower=30)
        flat = lc.flatten(window_length=481)
        t, y = flat.time.value, np.asarray(flat.flux.value, dtype=float)

        model = transitleastsquares(t, y)
        tls = model.power(
            period_min=max(0.52, r.period * 0.8),
            period_max=min(13.0, r.period * 1.2),
            R_star=float(tinfo["rad"]), M_star=float(tinfo["mass"]),
            show_progress_bar=False,
        )
        print(f"TIC {tic}: BLS P={r.period:.4f} | TLS P={tls.period:.4f} "
              f"SDE={tls.SDE:.1f} FAP={tls.FAP:.2e} snr={tls.snr:.1f} "
              f"depth={(1-tls.depth)*1e6:.0f} ppm odd/even={tls.odd_even_mismatch:.2f}sig "
              f"transits={tls.distinct_transit_count}")
        rows.append(dict(
            tic=tic, Tmag=r.Tmag, Teff=r.Teff, rad=r.rad,
            bls_period=r.period, bls_sde=r.sde, depth_ppm=r.depth_ppm,
            rp_earth=r.rp_earth, n_transits=r.n_transits,
            tls_period=round(float(tls.period), 5), tls_sde=round(float(tls.SDE), 2),
            tls_fap=float(tls.FAP), tls_snr=round(float(tls.snr), 2),
            tls_depth_ppm=round((1 - float(tls.depth)) * 1e6, 1),
            tls_oddeven_sigma=round(float(tls.odd_even_mismatch), 2),
            tls_transits=int(tls.distinct_transit_count),
        ))

    vetted = pd.DataFrame(rows)
    if len(vetted):
        vetted["passes_tls"] = (vetted["tls_sde"] >= 9) & (vetted["tls_fap"] < 1e-3)
        vetted = vetted.sort_values("tls_sde", ascending=False)
    vetted.to_csv(CAND / "vetted.csv", index=False)
    print(f"\nwrote {CAND / 'vetted.csv'}")
    if len(vetted):
        print(vetted[["tic", "bls_period", "tls_period", "tls_sde", "tls_fap",
                      "tls_snr", "rp_earth", "passes_tls"]].to_string(index=False))


if __name__ == "__main__":
    main()

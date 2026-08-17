"""Transit search over the unflagged sector-104 M/K-dwarf sample.

Per target: download the SPOC light curve, clean + flatten, BLS period
search, vet (transit count, odd/even depth, secondary eclipse), and record
everything. Candidates that survive get a 4-panel diagnostic figure.

Outputs:
  data/results.csv          one row per star searched (full sweep record)
  candidates/candidates.csv surviving candidates only
  figures/TIC<id>.png       diagnostics for each candidate
Resumable: already-searched TICs are skipped on rerun.
Run:  .venv/bin/python hunt.py
"""

import csv
import os
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import astropy.units as u
from astropy.timeseries import BoxLeastSquares
import lightkurve as lk

warnings.filterwarnings("ignore")

ROOT = Path(__file__).parent
DATA, FIGS, CAND = ROOT / "data", ROOT / "figures", ROOT / "candidates"
RESULTS = DATA / "results.csv"
CANDFILE = CAND / "candidates.csv"
CACHE = DATA / "lc_cache"
CACHE.mkdir(exist_ok=True)

SDE_MIN = 8.0        # BLS detection strength threshold
MIN_TRANSITS = 2
MAX_DEPTH = 0.10     # >10% dip = almost certainly an eclipsing binary
SUN_R_EARTH = 109.2

RES_COLS = ["tic", "Tmag", "Teff", "rad", "sde", "period", "t0", "depth_ppm",
            "duration_hr", "n_transits", "rp_earth", "oddeven_sigma",
            "secondary_sigma", "verdict"]


def in_transit_mask(time, period, t0, duration):
    phase = ((time - t0 + 0.5 * period) % period) - 0.5 * period
    return np.abs(phase) < 0.55 * duration


def analyze_one(row):
    tic = int(row.tic)
    fpath = CACHE / row.fname
    if not fpath.exists():
        urllib.request.urlretrieve(row.url, fpath)
    lc = lk.read(fpath, flux_column="pdcsap_flux")
    lc = lc.remove_nans().normalize().remove_outliers(sigma_upper=3, sigma_lower=30)
    flat = lc.flatten(window_length=481)
    binned = flat.bin(time_bin_size=10 * u.minute).remove_nans()
    t, y = binned.time.value, binned.flux.value
    if len(t) < 500:
        return dict(tic=tic, verdict="too_few_points")

    periods = np.linspace(0.5, 13.0, 40000)
    durations = np.array([0.03, 0.06, 0.10, 0.15])
    bls = BoxLeastSquares(t, y)
    res = bls.power(periods, durations)
    i = int(np.argmax(res.power))
    period, t0 = float(res.period[i]), float(res.transit_time[i])
    depth, duration = float(res.depth[i]), float(res.duration[i])
    sde = float((res.power.max() - np.median(res.power)) / np.std(res.power))

    itr = in_transit_mask(t, period, t0, duration)
    oot = ~in_transit_mask(t, period, t0, 2.5 * duration)
    scatter = np.std(y[oot])

    # distinct transit epochs actually covered by data
    epochs = np.unique(np.round((t[itr] - t0) / period).astype(int))
    n_transits = len(epochs)

    # odd/even depth consistency (eclipsing binaries differ; planets don't)
    ep_idx = np.round((t - t0) / period).astype(int)
    d_odd = 1 - np.mean(y[itr & (ep_idx % 2 == 1)]) if np.any(itr & (ep_idx % 2 == 1)) else np.nan
    d_even = 1 - np.mean(y[itr & (ep_idx % 2 == 0)]) if np.any(itr & (ep_idx % 2 == 0)) else np.nan
    n_odd = max(np.sum(itr & (ep_idx % 2 == 1)), 1)
    n_even = max(np.sum(itr & (ep_idx % 2 == 0)), 1)
    oe_err = scatter * np.sqrt(1 / n_odd + 1 / n_even)
    oe_sigma = float(abs(d_odd - d_even) / oe_err) if np.isfinite(d_odd + d_even) else 0.0

    # secondary eclipse at phase 0.5 (planet: none; binary: often yes)
    sec = in_transit_mask(t, period, t0 + 0.5 * period, duration)
    sec_depth = 1 - np.mean(y[sec]) if np.any(sec) else 0.0
    sec_sigma = float(sec_depth / (scatter / np.sqrt(max(np.sum(sec), 1))))

    rp = float(np.sqrt(max(depth, 0)) * row.rad * SUN_R_EARTH) if row.rad > 0 else np.nan

    verdict = "no_signal"
    if sde >= SDE_MIN and depth > 0 and n_transits >= MIN_TRANSITS:
        if depth > MAX_DEPTH:
            verdict = "eclipsing_binary_depth"
        elif oe_sigma > 5:
            verdict = "eclipsing_binary_oddeven"
        elif sec_sigma > 5:
            verdict = "eclipsing_binary_secondary"
        else:
            verdict = "CANDIDATE"

    out = dict(tic=tic, Tmag=row.Tmag, Teff=row.Teff, rad=row.rad, sde=round(sde, 2),
               period=round(period, 5), t0=round(t0, 5), depth_ppm=round(depth * 1e6, 1),
               duration_hr=round(duration * 24, 2), n_transits=n_transits,
               rp_earth=round(rp, 2) if np.isfinite(rp) else "",
               oddeven_sigma=round(oe_sigma, 2), secondary_sigma=round(sec_sigma, 2),
               verdict=verdict)

    if verdict == "CANDIDATE":
        folded = flat.fold(period=period * u.day, epoch_time=t0)
        fig, axes = plt.subplots(2, 2, figsize=(13, 8))
        lc.scatter(ax=axes[0, 0], s=1)
        axes[0, 0].set_title(f"TIC {tic} — raw (Tmag={row.Tmag:.1f}, R*={row.rad:.2f} Rsun)")
        flat.scatter(ax=axes[0, 1], s=1)
        axes[0, 1].set_title("flattened")
        axes[1, 0].plot(res.period, res.power, lw=0.4)
        axes[1, 0].axvline(period, color="r", ls="--", alpha=0.6)
        axes[1, 0].set_title(f"BLS — P={period:.4f} d, SDE={sde:.1f}")
        axes[1, 0].set_xlabel("period [d]")
        folded.scatter(ax=axes[1, 1], s=1)
        fb = folded.bin(time_bin_size=10 * u.minute)
        axes[1, 1].plot(fb.time.value, fb.flux.value, color="r", lw=1.2)
        axes[1, 1].set_xlim(-4 * duration, 4 * duration)
        axes[1, 1].set_title(f"folded — depth {depth*1e6:.0f} ppm, Rp≈{rp:.1f} R_E, {n_transits} transits")
        fig.tight_layout()
        fig.savefig(FIGS / f"TIC{tic}.png", dpi=130)
        plt.close(fig)

    fpath.unlink(missing_ok=True)  # keep the cache small
    return out


def append_row(path, row, cols):
    new = not path.exists()
    with open(path, "a", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        if new:
            w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    targets = pd.read_csv(DATA / "targets.csv")
    done = set()
    if RESULTS.exists():
        done = set(pd.read_csv(RESULTS)["tic"].astype(int))
    todo = targets[~targets["tic"].isin(done)]
    print(f"{len(todo)} targets to search ({len(done)} already done)")
    n_cand = 0
    for k, row in enumerate(todo.itertuples(), 1):
        try:
            out = analyze_one(row)
        except Exception as e:
            out = dict(tic=int(row.tic), verdict=f"error:{type(e).__name__}")
        for c in RES_COLS:
            out.setdefault(c, "")
        append_row(RESULTS, out, RES_COLS)
        if out["verdict"] == "CANDIDATE":
            n_cand += 1
            append_row(CANDFILE, out, RES_COLS)
            print(f"  *** CANDIDATE: TIC {out['tic']}  P={out['period']} d  "
                  f"depth={out['depth_ppm']} ppm  SDE={out['sde']}")
        if k % 20 == 0:
            print(f"  [{k}/{len(todo)}] searched, {n_cand} candidates so far")
    print(f"done: {len(todo)} searched, {n_cand} candidates")

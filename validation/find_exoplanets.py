"""Find exoplanets in real TESS data with a Box Least Squares transit search.

For each target star:
  1. Download the 2-min cadence SPOC light curve from the MAST archive
  2. Clean (remove NaNs / outliers) and flatten (remove stellar variability)
  3. Run a BLS periodogram to search for periodic transit-shaped dips
  4. Report period, epoch, depth, duration, S/N and an estimated planet radius
  5. Save a 4-panel diagnostic figure per target in figures/

Run:  .venv/bin/python find_exoplanets.py
"""

import warnings
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightkurve as lk
import astropy.units as u

warnings.filterwarnings("ignore")

FIGDIR = Path(__file__).parent / "figures"
FIGDIR.mkdir(exist_ok=True)

# name, max sectors to use, known period (days) for validation, notes
TARGETS = [
    ("WASP-18", 1, 0.9415, "hot Jupiter, ~1% deep transit"),
    ("LHS 3844", 1, 0.4629, "ultra-short-period rocky planet"),
    ("Pi Mensae", 2, 6.2679, "shallow ~300 ppm mini-Neptune"),
]

SUN_R_EARTH = 109.2  # Earth radii per solar radius


def analyze(name, n_sectors, known_period, note):
    print(f"\n=== {name} ({note}) ===")
    search = lk.search_lightcurve(name, mission="TESS", author="SPOC", exptime=120)
    if len(search) == 0:
        print("  no SPOC light curve found")
        return None
    lcs = search[:n_sectors].download_all()
    lc = lcs.stitch().remove_nans().remove_outliers(sigma_upper=3, sigma_lower=20)
    star_radius = lc.meta.get("RADIUS")  # solar radii, from the TESS Input Catalog

    # flatten with a window much longer than a transit so dips survive
    flat = lc.flatten(window_length=401)

    # astropy BLS directly (lightkurve's wrapper has a unit bug with explicit grids)
    from astropy.timeseries import BoxLeastSquares

    t, y = flat.time.value, flat.flux.value
    period_grid = np.linspace(0.3, 15, 60000)
    bls_res = BoxLeastSquares(t, y).power(period_grid, np.array([0.05, 0.1, 0.2]))
    best = int(np.argmax(bls_res.power))
    period = bls_res.period[best] * u.day
    t0 = bls_res.transit_time[best]
    depth = float(bls_res.depth[best])
    duration = bls_res.duration[best] * u.day

    # detection strength: peak height vs noise floor of the BLS power spectrum
    power = bls_res.power
    sde = (power.max() - np.median(power)) / np.std(power)

    rp_earth = None
    if star_radius and depth > 0:
        rp_earth = float(np.sqrt(depth) * star_radius * SUN_R_EARTH)

    folded = flat.fold(period=period, epoch_time=t0)

    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    lc.scatter(ax=axes[0, 0], s=1)
    axes[0, 0].set_title(f"{name} — raw TESS light curve")
    flat.scatter(ax=axes[0, 1], s=1)
    axes[0, 1].set_title("flattened")
    axes[1, 0].plot(bls_res.period, power, lw=0.5)
    axes[1, 0].axvline(period.value, color="r", ls="--", alpha=0.6)
    axes[1, 0].set_xlabel("period [d]")
    axes[1, 0].set_ylabel("BLS power")
    axes[1, 0].set_title(f"BLS periodogram — peak at {period.value:.4f} d")
    folded.scatter(ax=axes[1, 1], s=1)
    binned = folded.bin(time_bin_size=5 * u.minute)
    axes[1, 1].plot(binned.time.value, binned.flux.value, color="r", lw=1.2)
    axes[1, 1].set_xlim(-0.3, 0.3)
    axes[1, 1].set_title("phase-folded transit (red = 5-min bins)")
    fig.tight_layout()
    outfile = FIGDIR / f"{name.replace(' ', '_')}.png"
    fig.savefig(outfile, dpi=130)
    plt.close(fig)

    result = dict(
        name=name,
        period=float(period.value),
        known_period=known_period,
        t0=float(t0),
        depth_ppm=float(depth) * 1e6,
        duration_hr=float(duration.to(u.hour).value),
        sde=float(sde),
        rp_earth=rp_earth,
        figure=str(outfile),
    )
    match = abs(result["period"] - known_period) / known_period < 0.01
    print(f"  detected period : {result['period']:.4f} d  (known: {known_period} d)"
          f"  {'MATCH' if match else 'MISMATCH'}")
    print(f"  transit depth   : {result['depth_ppm']:.0f} ppm")
    print(f"  duration        : {result['duration_hr']:.2f} h")
    print(f"  detection S/N   : {result['sde']:.1f}")
    if rp_earth:
        print(f"  planet radius   : {rp_earth:.2f} R_Earth  (star: {star_radius:.2f} R_Sun)")
    print(f"  figure          : {outfile}")
    return result


if __name__ == "__main__":
    results = [analyze(*t) for t in TARGETS]
    found = [r for r in results if r and r["sde"] > 7]
    print(f"\n{len(found)}/{len(TARGETS)} targets show a significant transit signal (S/N > 7).")

"""Multi-sector follow-up of the two marginal candidates: stitch every
available SPOC sector and re-run TLS near the candidate period."""
import warnings
warnings.filterwarnings("ignore")
import numpy as np
import lightkurve as lk
from transitleastsquares import transitleastsquares

CANDS = [(260817968, 1.261), (87329149, 3.336)]

def main():
    for tic, p0 in CANDS:
        sr = lk.search_lightcurve(f"TIC {tic}", mission="TESS", author="SPOC", exptime=120)
        print(f"\nTIC {tic}: {len(sr)} SPOC sectors: {list(sr.mission)}")
        if len(sr) < 2:
            print("  no archival sectors — cannot strengthen; needs future data")
            continue
        lc = sr.download_all().stitch().remove_nans().remove_outliers(sigma_upper=3, sigma_lower=30)
        flat = lc.flatten(window_length=481)
        t, y = flat.time.value, np.asarray(flat.flux.value, dtype=float)
        tls = transitleastsquares(t, y).power(
            period_min=p0 * 0.95, period_max=p0 * 1.05, show_progress_bar=False)
        print(f"  multi-sector TLS: P={tls.period:.5f} SDE={tls.SDE:.1f} "
              f"FAP={tls.FAP:.2e} snr={tls.snr:.1f} depth={(1-tls.depth)*1e6:.0f} ppm "
              f"transits={tls.distinct_transit_count}")

if __name__ == "__main__":
    main()

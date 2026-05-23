"""
run_bayesmsd.py
---------------
Load MATLAB tracks (cell array saved as .mat) and fit fBM parameters with
bayesmsd, then compare to MATLAB fitFBM_MLE results.

Usage
-----
python run_bayesmsd.py tracks.mat \
    --dt 0.02 \
    --exposure_fraction 0.5 \
    --n_tracks 500 \
    --min_length 15 \
    --matlab_K 0.0032 \
    --matlab_alpha 0.48 \
    --matlab_sigma 0.025

The .mat file should contain a variable named 'tracks' that is a MATLAB
cell array (N x 1 or 1 x N), where each cell is an (L x 2) matrix of
[x, y] positions in µm.

Convention note
---------------
bayesmsd fits MSD_1d(τ) = Γ · τ^α per dimension.
MATLAB fitFBM_MLE fits with Var[Δx]_1d = 2·K·τ^α (no-blur limit).
Therefore:  Γ_bayesmsd = 2 · K_MATLAB
Localization noise:  bayesmsd σ² (per dim) = sigma_MATLAB²
"""

import argparse
import sys
import warnings
import numpy as np

# ── MATLAB file loading ───────────────────────────────────────────────────────

def load_mat_tracks(path, var_name='tracks'):
    """
    Load a MATLAB cell array of tracks from a .mat file.

    Tries scipy.io.loadmat (v5/v6) first; falls back to h5py (v7.3 / -v7.3).

    Returns a list of numpy arrays, each shape (L, 2) with columns [x, y] in µm.
    """
    tracks = None

    # --- Try scipy (MATLAB v5/v6) ---
    try:
        import scipy.io
        mat = scipy.io.loadmat(path, squeeze_me=True, struct_as_record=False)
        if var_name not in mat:
            raise KeyError(f"Variable '{var_name}' not found. "
                           f"Available: {[k for k in mat if not k.startswith('_')]}")
        raw = mat[var_name]
        # squeeze_me flattens 1xN or Nx1 cell arrays to a 1-D object array
        if raw.dtype == object:
            tracks = [np.asarray(raw[i], dtype=float) for i in range(raw.size)]
        else:
            raise ValueError("Expected a MATLAB cell array (dtype=object).")
        return tracks
    except NotImplementedError:
        pass  # v7.3 – fall through to h5py

    # --- Try h5py (MATLAB v7.3 / HDF5) ---
    try:
        import h5py
    except ImportError:
        sys.exit("Install h5py to read MATLAB v7.3 files:  pip install h5py")

    tracks = []
    with h5py.File(path, 'r') as f:
        if var_name not in f:
            raise KeyError(f"Variable '{var_name}' not found in HDF5 file.")
        cell = f[var_name]
        # HDF5 cell arrays are stored as object arrays of references
        refs = cell[()]
        refs = refs.flatten()
        for ref in refs:
            arr = f[ref][()].T  # HDF5 stores column-major; transpose to (L, d)
            tracks.append(np.asarray(arr, dtype=float))

    return tracks


# ── Preprocessing ─────────────────────────────────────────────────────────────

def filter_tracks(tracks, min_length=15, min_step_var=0.0):
    """Remove short and near-immobile tracks."""
    kept = []
    for tr in tracks:
        if tr.ndim != 2 or tr.shape[1] < 2:
            continue
        if tr.shape[0] < min_length:
            continue
        if min_step_var > 0:
            step_var = np.mean(np.diff(tr[:, 0])**2 + np.diff(tr[:, 1])**2)
            if step_var < min_step_var:
                continue
        kept.append(tr[:, :2])   # keep only x, y columns
    return kept


def subsample(tracks, n, seed=42):
    """Randomly subsample n tracks (without replacement)."""
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(tracks), size=min(n, len(tracks)), replace=False)
    return [tracks[i] for i in sorted(idx)]


# ── bayesmsd fit ──────────────────────────────────────────────────────────────

def run_bayesmsd(tracks, dt, exposure_fraction, ci=True, verbosity=1):
    """
    Fit fBM MSD (powerlaw + localization noise) to a list of 2D tracks.

    Parameters
    ----------
    tracks           : list of (L, 2) arrays  [x, y] in µm
    dt               : frame interval in seconds
    exposure_fraction: Te / dt  (0 = stroboscopic, 1 = full-frame)
    ci               : compute profile-likelihood 95% CIs
    verbosity        : 0 = silent, 1 = progress

    Returns
    -------
    result  : dict  raw bayesmsd result
    mci     : dict  {param: (estimate, [lo, hi])}  or None if ci=False
    fit     : NPXFit object
    """
    try:
        from bayesmsd.lib import NPXFit
        from bayesmsd import Profiler
    except ImportError as e:
        sys.exit(f"bayesmsd import error: {e}\nInstall with:  pip install bayesmsd")

    # bayesmsd works in units of the frame interval internally.
    # Scale tracks from µm to the same units (µm kept; dt in seconds passed
    # via the data lag structure — bayesmsd treats consecutive rows as lag-1).
    # There is no explicit dt argument; lag times are integers × dt when you
    # query the MSD later.

    # motion_blur_f = Te / dt = exposure_fraction
    # parametrization '(log(αΓ), α)' matches our log(K·α) parametrization and
    # reduces K-α posterior correlation (same reasoning as in MATLAB).
    fit = NPXFit(
        tracks,
        ss_order=1,             # non-stationary increments — correct for fBM
        n=0,                    # pure powerlaw (no spline extension)
        motion_blur_f=exposure_fraction,
        parametrization='(log(αΓ), α)',
    )

    if verbosity:
        print(f"Running bayesmsd NPXFit on {len(tracks)} tracks ...")

    result = fit.run(verbosity=verbosity)

    mci = None
    if ci:
        if verbosity:
            print("Computing profile-likelihood CIs ...")
        profiler = Profiler(fit, profiling=True, conf=0.95)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mci = profiler.find_MCI()

    return result, mci, fit


# ── Parameter extraction and unit conversion ──────────────────────────────────

def extract_params(result, mci, dt):
    """
    Pull fitted parameters and convert to MATLAB convention.

    bayesmsd units (lag = integer frames; MSD in µm²):
      Γ  = prefactor of 1D MSD:  MSD_1d(τ) = Γ · τ^α  [µm² / frame^α]
      α  = anomalous exponent
      σ² = localization noise variance per dimension [µm²]

    MATLAB units:
      K  = generalized diffusion coefficient  [µm² / s^α]
         = Γ / (2 · dt^α)      (factor 2 from Var[Δx] = 2K·τ^α convention;
                                  dt^α converts frame^α → s^α)
      α  = same
      σ  = localization precision [µm] = sqrt(σ²_bayesmsd)
    """
    params = result['params']
    out = {}

    # α
    alpha = params['α']
    out['alpha'] = alpha

    # Γ from log(αΓ) parametrization: αΓ = exp(log(αΓ)), so Γ = αΓ / α
    # For dim 0:
    logaG_key = 'log(αΓ) (dim 0)'
    if logaG_key in params:
        aG = np.exp(params[logaG_key])
        G  = aG / alpha
    else:
        # fallback: log(Γ) parametrization
        G = np.exp(params.get('log(Γ) (dim 0)', np.nan))

    # Convert Γ [µm²/frame^α] → K [µm²/s^α]
    # MSD_1d = Γ · n^α  (n in frames)  = Γ/dt^α · (n·dt)^α  = (Γ/dt^α) · τ^α
    # MATLAB: MSD_1d = 2K·τ^α  →  K = Γ / (2 · dt^α)
    out['K']  = G / (2.0 * dt**alpha)
    out['Ka'] = out['K'] * alpha

    # σ per dimension (bayesmsd σ² is variance per dim)
    sig2_key = 'log(σ²) (dim 0)'
    if sig2_key in params:
        out['sigma'] = np.sqrt(np.exp(params[sig2_key]))
    else:
        out['sigma'] = np.nan

    out['loglik'] = result['logL']

    # CIs
    if mci is not None:
        def _ci(key):
            if key in mci:
                _, (lo, hi) = mci[key]
                return lo, hi
            return np.nan, np.nan

        a_lo, a_hi     = _ci('α')
        out['alpha_CI'] = (a_lo, a_hi)

        # Propagate CI on log(αΓ) → K  (approximate, ignoring α–Γ correlation)
        logaG_lo, logaG_hi = _ci(logaG_key)
        if not np.isnan(logaG_lo):
            K_lo = np.exp(logaG_lo) / alpha / (2.0 * dt**alpha)
            K_hi = np.exp(logaG_hi) / alpha / (2.0 * dt**alpha)
            out['K_CI'] = (K_lo, K_hi)
        else:
            out['K_CI'] = (np.nan, np.nan)

        sig2_lo, sig2_hi = _ci(sig2_key)
        out['sigma_CI'] = (np.sqrt(np.exp(sig2_lo)), np.sqrt(np.exp(sig2_hi))) \
                          if not np.isnan(sig2_lo) else (np.nan, np.nan)

    return out


# ── Reporting ─────────────────────────────────────────────────────────────────

def print_comparison(bms, matlab_K=None, matlab_alpha=None, matlab_sigma=None):
    """Print side-by-side comparison table."""
    sep = '─' * 62

    print()
    print(sep)
    print(f"{'':20s}  {'bayesmsd':>18s}  {'MATLAB MLE':>18s}")
    print(sep)

    def row(name, bval, mval, ci=None):
        bstr = f"{bval:.4f}"
        if ci:
            bstr += f"  [{ci[0]:.4f}, {ci[1]:.4f}]"
        mstr = f"{mval:.4f}" if mval is not None else "—"
        print(f"  {name:<18s}  {bstr:>18s}  {mstr:>18s}")

    row("alpha",  bms['alpha'], matlab_alpha,
        bms.get('alpha_CI'))
    row("K  (µm²/s^α)", bms['K'], matlab_K,
        bms.get('K_CI'))
    row("Ka = K·α", bms['Ka'], (matlab_K * matlab_alpha) if matlab_K and matlab_alpha else None)
    row("sigma (µm)", bms['sigma'], matlab_sigma,
        bms.get('sigma_CI'))

    print(sep)
    print(f"  {'log-likelihood':<18s}  {bms['loglik']:>18.2f}")
    print(sep)

    print()
    print("Convention note:")
    print("  bayesmsd MSD_1d(τ) = Γ·τ^α  vs  MATLAB Var[Δx]_1d = 2K·τ^α")
    print("  K_MATLAB = Γ / (2·dt^α);  σ_MATLAB = sqrt(σ²_bayesmsd per dim)")
    print()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Compare bayesmsd fBM fit to MATLAB fitFBM_MLE on the same tracks.")
    p.add_argument("mat_file",
                   help="Path to .mat file containing the 'tracks' cell array.")
    p.add_argument("--var",            default="tracks",
                   help="MATLAB variable name for the cell array (default: tracks).")
    p.add_argument("--dt",             type=float, required=True,
                   help="Frame interval in seconds.")
    p.add_argument("--exposure_fraction", type=float, default=1.0,
                   help="Te / dt exposure fraction (default: 1.0).")
    p.add_argument("--n_tracks",       type=int,   default=500,
                   help="Number of tracks to subsample (default: 500).")
    p.add_argument("--min_length",     type=int,   default=15,
                   help="Minimum track length in frames (default: 15).")
    p.add_argument("--min_step_var",   type=float, default=0.0,
                   help="Minimum mean squared step size µm² to include track "
                        "(default: 0 = no filter).  Use ~4*sigma^2 to drop "
                        "stuck particles.")
    p.add_argument("--seed",           type=int,   default=42,
                   help="Random seed for subsampling (default: 42).")
    p.add_argument("--no_ci",          action="store_true",
                   help="Skip profile-likelihood CIs (faster).")
    # Optional MATLAB results for comparison
    p.add_argument("--matlab_K",       type=float, default=None,
                   help="MATLAB fitFBM_MLE K estimate (µm²/s^α) for comparison.")
    p.add_argument("--matlab_alpha",   type=float, default=None,
                   help="MATLAB fitFBM_MLE alpha estimate for comparison.")
    p.add_argument("--matlab_sigma",   type=float, default=None,
                   help="MATLAB fitFBM_MLE sigma estimate (µm) for comparison.")
    args = p.parse_args()

    # --- Load ---
    print(f"Loading '{args.var}' from {args.mat_file} ...")
    try:
        tracks = load_mat_tracks(args.mat_file, var_name=args.var)
    except Exception as e:
        sys.exit(f"Error loading file: {e}")
    print(f"  Loaded {len(tracks)} tracks.")

    # --- Filter ---
    tracks = filter_tracks(tracks,
                           min_length=args.min_length,
                           min_step_var=args.min_step_var)
    print(f"  {len(tracks)} tracks after length/step-var filter.")

    if len(tracks) == 0:
        sys.exit("No tracks remain after filtering.")

    # --- Subsample ---
    if args.n_tracks < len(tracks):
        tracks = subsample(tracks, args.n_tracks, seed=args.seed)
        print(f"  Subsampled to {len(tracks)} tracks (seed={args.seed}).")

    lengths = [t.shape[0] for t in tracks]
    print(f"  Track lengths: min={min(lengths)}, median={int(np.median(lengths))}, "
          f"max={max(lengths)}")

    # --- Fit ---
    result, mci, fit = run_bayesmsd(
        tracks,
        dt=args.dt,
        exposure_fraction=args.exposure_fraction,
        ci=not args.no_ci,
    )

    # --- Extract + convert ---
    bms = extract_params(result, mci, args.dt)

    # --- Report ---
    print_comparison(bms,
                     matlab_K=args.matlab_K,
                     matlab_alpha=args.matlab_alpha,
                     matlab_sigma=args.matlab_sigma)

    # --- Raw bayesmsd params for reference ---
    print("Raw bayesmsd parameters:")
    for k, v in result['params'].items():
        ci_str = ""
        if mci and k in mci:
            _, (lo, hi) = mci[k]
            ci_str = f"  95% CI: [{lo:.4f}, {hi:.4f}]"
        print(f"  {k:30s} = {v:.4f}{ci_str}")
    print()


if __name__ == "__main__":
    main()

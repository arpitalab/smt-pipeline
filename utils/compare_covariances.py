"""
compare_covariances.py
----------------------
Compare the fBM displacement covariance matrix from:
  1. MATLAB fitFBM_MLE  (Backlund et al. exact Toeplitz kernel)
  2. bayesmsd NPXFit    (GP with ss_order=1 + imaging() decorator)

for identical parameters (K, alpha, sigma, dt, exposure_fraction).

The script:
  - Builds both covariance matrices row-by-row (all lags 0..nd-1)
  - Shows where and by how much they differ
  - Identifies whether the discrepancy can explain the systematic
    ~0.11 offset in fitted alpha between the two methods.

Usage
-----
python compare_covariances.py \
    --K 0.008 --alpha 0.57 --sigma 0.054 \
    --dt 0.2 --exposure_fraction 0.05 --nd 19

All parameters match the (per-subtrack) convention used in:
  MATLAB:   Var[Δx]_1d = 2K τ^α  plus noise sigma
  bayesmsd: MSD_1d(τ) = Γ τ^α   plus noise σ² per dim
  Relation: K = Γ / (2·dt^α),  sigma = sqrt(noise2)
"""

import argparse
import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# 1.  MATLAB Backlund exact kernel
# ─────────────────────────────────────────────────────────────────────────────

def psi_fbm_matlab(tau, Te, alpha):
    """
    Exposure-time-averaged fBM position covariance kernel (Backlund et al. 2015).

      Psi(tau, Te, alpha) = [(tau+Te)^(a+2) + |tau-Te|^(a+2) - 2*tau^(a+2)]
                            / [Te^2 * (a+1) * (a+2)]

    Limit Te→0 (stroboscopic): Psi → tau^alpha.
    """
    if Te < 1e-12:
        return tau ** alpha
    return (
        (tau + Te) ** (alpha + 2) + abs(tau - Te) ** (alpha + 2) - 2 * tau ** (alpha + 2)
    ) / (Te ** 2 * (alpha + 1) * (alpha + 2))


def build_cov_matlab(nd, K, alpha, sigma, dt, Te):
    """
    Build the (nd x nd) Toeplitz displacement covariance matrix.

    The n-th off-diagonal (0-indexed) is:
      c[n] = K * [Psi((n+1)*dt) + Psi(|n-1|*dt) - 2*Psi(n*dt)]
    with
      c[0] += 2*sigma^2    (diagonal noise)
      c[1] -= sigma^2      (nearest-neighbour anti-correlation)
    """
    c = np.zeros(nd)
    for n in range(nd):
        psi_p = psi_fbm_matlab((n + 1) * dt, Te, alpha)
        psi_0 = psi_fbm_matlab(n * dt, Te, alpha)  # Psi(0) = 2Te^alpha/((a+1)(a+2))
        psi_m = psi_fbm_matlab(abs(n - 1) * dt, Te, alpha)
        c[n] = K * (psi_p + psi_m - 2 * psi_0)

    c[0] += 2 * sigma ** 2
    if nd > 1:
        c[1] -= sigma ** 2

    from scipy.linalg import toeplitz
    return toeplitz(c), c


# ─────────────────────────────────────────────────────────────────────────────
# 2.  bayesmsd covariance (ss_order=1, imaging decorator)
# ─────────────────────────────────────────────────────────────────────────────

def msd_imaging(n_frames, G, alpha, f, noise2):
    """
    bayesmsd MSD_eff for a powerlaw fBM with imaging correction.

    Works in frame units (lag = integer frames, G in µm²/frame^alpha).

      MSD_raw(n) = G * n^alpha
      B          = MSD_raw(f) / ((alpha+1)*(alpha+2))   [f = Te/dt_frame]
      phi        = f / n
      b(n)       = ((1+phi)^(a+2) + (1-phi)^(a+2) - 2) / (phi^2*(a+1)*(a+2))
      MSD_eff(n) = b(n) * MSD_raw(n) - 2*B + 2*noise2

    The MSDfun decorator in bayesmsd forces MSD_eff(0) = 0 (handled separately).

    Parameters
    ----------
    n_frames : float or array  (lag in frames, may include 0)
    G        : float           prefactor [µm²/frame^alpha]
    alpha    : float
    f        : float           Te/dt_frame (exposure fraction, 0=stroboscopic)
    noise2   : float           localization noise variance per dim [µm²]

    Returns
    -------
    msd_eff  : array
    """
    n = np.atleast_1d(np.asarray(n_frames, dtype=float))
    out = np.zeros_like(n)

    if f == 0:
        # Stroboscopic: no motion blur
        pos = n > 0
        out[pos] = G * n[pos] ** alpha + 2 * noise2
        return out

    B = G * f ** alpha / ((alpha + 1) * (alpha + 2))

    pos = n > 0
    phi = f / n[pos]
    b = ((1 + phi) ** (alpha + 2) + (1 - phi) ** (alpha + 2) - 2) / (
        phi ** 2 * (alpha + 1) * (alpha + 2)
    )
    out[pos] = b * G * n[pos] ** alpha - 2 * B + 2 * noise2
    # n==0: MSDfun forces 0 (already initialised to 0)

    return out


def build_cov_bayesmsd(nd, G, alpha, f, noise2):
    """
    Build the (nd x nd) increment covariance for ss_order=1.

    The gp.msd2C_fun formula (for integer lags i, j with |i-j| = n):
      C[i,j] = 0.5 * (MSD(n+1) + MSD(|n-1|) - 2*MSD(n))

    where MSD is MSD_eff (including motion blur + noise) evaluated at
    integer frame lags, and MSD(0) = 0 by convention.
    """
    # Precompute MSD_eff at lags 0 .. nd+1
    lags = np.arange(nd + 2, dtype=float)
    msd = msd_imaging(lags, G, alpha, f, noise2)
    msd[0] = 0.0   # enforced by MSDfun decorator

    c = np.zeros(nd)
    for n in range(nd):
        c[n] = 0.5 * (msd[n + 1] + msd[abs(n - 1)] - 2 * msd[n])

    from scipy.linalg import toeplitz
    return toeplitz(c), c


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Analytic cross-check: are the first rows equal?
# ─────────────────────────────────────────────────────────────────────────────

def first_row_comparison(K, alpha, sigma, dt, frac, nd):
    """
    Compare the first row (= all unique lags) of both covariance matrices.
    Returns (c_matlab, c_bayesmsd, rel_diff).
    """
    Te = frac * dt

    # bayesmsd convention: K = G / (2 * dt^alpha)  →  G = 2K * dt^alpha
    G      = 2.0 * K * dt ** alpha
    noise2 = sigma ** 2
    f      = frac  # = Te / dt_frame

    _, c_mat = build_cov_matlab(nd, K, alpha, sigma, dt, Te)
    _, c_bms = build_cov_bayesmsd(nd, G, alpha, f, noise2)

    return c_mat, c_bms


# ─────────────────────────────────────────────────────────────────────────────
# 4.  What fitted alpha would bayesmsd recover if MATLAB's true matrix is used?
# ─────────────────────────────────────────────────────────────────────────────

def fit_alpha_from_cov(c_target, K0, alpha0, sigma0, dt, frac, nd,
                       method='bms', n_grid=200, tol=1e-6):
    """
    Given a target first-row c_target (true covariance), find the alpha that
    minimises ||c_model(alpha) - c_target||² when K and sigma are also varied.

    This is a rough 1-D grid search over alpha (optimising K, sigma analytically
    is hard, so we scan a 3-D grid of alpha only and do joint minimisation via
    scipy).

    Returns fitted (K, alpha, sigma) and the residual norm.
    """
    from scipy.optimize import minimize

    Te = frac * dt

    def cost(params):
        K_, logit_a, log_s = params
        a_    = 2 / (1 + np.exp(-logit_a))
        s_    = np.exp(log_s)
        K_    = np.exp(K_)
        if method == 'bms':
            G_   = 2 * K_ * dt ** a_
            f_   = frac
            _, c = build_cov_bayesmsd(nd, G_, a_, f_, s_ ** 2)
        else:
            _, c = build_cov_matlab(nd, K_, a_, s_, dt, Te)
        return np.sum((c - c_target) ** 2)

    th0 = [np.log(K0), np.log((alpha0 / 2) / (1 - alpha0 / 2)), np.log(sigma0)]
    res = minimize(cost, th0, method='Nelder-Mead',
                   options={'xatol': tol, 'fatol': tol ** 2, 'maxiter': 50000})

    K_fit    = np.exp(res.x[0])
    alpha_fit = 2 / (1 + np.exp(-res.x[1]))
    sigma_fit = np.exp(res.x[2])
    return K_fit, alpha_fit, sigma_fit, res.fun


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="Compare MATLAB vs bayesmsd fBM covariance matrices.")
    p.add_argument("--K",    type=float, default=0.0079,
                   help="Diffusion coeff [µm²/s^alpha] (default from tracks_culled fit).")
    p.add_argument("--alpha", type=float, default=0.57,
                   help="Anomalous exponent (default: 0.57).")
    p.add_argument("--sigma", type=float, default=0.054,
                   help="Localization precision [µm] (default: 0.054).")
    p.add_argument("--dt",   type=float, default=0.2,
                   help="Frame interval [s] (default: 0.2).")
    p.add_argument("--exposure_fraction", type=float, default=0.05,
                   help="Te/dt (default: 0.05).")
    p.add_argument("--nd",   type=int,   default=19,
                   help="Number of displacements per subtrack = subtrack_length-1 (default: 19).")
    p.add_argument("--cross_fit", action="store_true",
                   help="Also fit alpha to the MATLAB covariance using the bayesmsd model "
                        "(shows how much alpha shifts purely from covariance model mismatch).")
    args = p.parse_args()

    K     = args.K
    alpha = args.alpha
    sigma = args.sigma
    dt    = args.dt
    frac  = args.exposure_fraction
    nd    = args.nd
    Te    = frac * dt
    G     = 2.0 * K * dt ** alpha   # bayesmsd Γ [µm²/frame^alpha]

    print()
    print("=" * 70)
    print(" fBM Covariance Matrix Comparison:  MATLAB  vs  bayesmsd")
    print("=" * 70)
    print(f"  Parameters:  K={K}  alpha={alpha}  sigma={sigma}")
    print(f"               dt={dt} s  frac={frac}  Te={Te:.4f} s  nd={nd}")
    print(f"  bayesmsd G = 2·K·dt^alpha = {G:.6f} µm²/frame^alpha")
    print()

    c_mat, c_bms = first_row_comparison(K, alpha, sigma, dt, frac, nd)

    # Relative difference (avoid divide-by-zero on tiny off-diagonals)
    eps    = 1e-14
    c_abs  = np.maximum(0.5 * (np.abs(c_mat) + np.abs(c_bms)), eps)
    reldiff = (c_bms - c_mat) / c_abs * 100  # %

    print(f"  {'lag (n)':<8}  {'MATLAB c(n)':>14}  {'bayesmsd c(n)':>14}  "
          f"{'diff':>14}  {'rel diff %':>10}")
    print("  " + "-" * 65)
    for n in range(nd):
        diff = c_bms[n] - c_mat[n]
        print(f"  {n:<8}  {c_mat[n]:>14.6e}  {c_bms[n]:>14.6e}  "
              f"{diff:>14.3e}  {reldiff[n]:>10.4f}")

    print()
    max_rel = np.max(np.abs(reldiff))
    rms_rel = np.sqrt(np.mean(reldiff ** 2))
    print(f"  Max |rel diff|: {max_rel:.4f}%   RMS rel diff: {rms_rel:.4f}%")

    # Check positive-definiteness of both matrices
    from scipy.linalg import toeplitz
    Sigma_mat = toeplitz(c_mat)
    Sigma_bms = toeplitz(c_bms)

    ev_mat = np.linalg.eigvalsh(Sigma_mat)
    ev_bms = np.linalg.eigvalsh(Sigma_bms)
    print(f"  MATLAB Sigma:  min eigenvalue = {ev_mat.min():.4e}  (PD: {ev_mat.min()>0})")
    print(f"  bayesmsd Sigma: min eigenvalue = {ev_bms.min():.4e}  (PD: {ev_bms.min()>0})")

    # ---- Stroboscopic (f=0) comparison ----------------------------------------
    print()
    print("─" * 70)
    print("  Stroboscopic limit (f=0): both should be identical")
    c_mat0, c_bms0 = first_row_comparison(K, alpha, sigma, dt, 0.0, nd)
    reldiff0 = (c_bms0 - c_mat0) / np.maximum(0.5 * (np.abs(c_mat0) + np.abs(c_bms0)), eps) * 100
    max_rel0 = np.max(np.abs(reldiff0))
    rms_rel0 = np.sqrt(np.mean(reldiff0 ** 2))
    print(f"  Max |rel diff|: {max_rel0:.6f}%   RMS rel diff: {rms_rel0:.6f}%")
    if max_rel0 < 1e-8:
        print("  → Stroboscopic covariances are numerically identical ✓")
    else:
        print("  → WARNING: stroboscopic covariances differ — check implementation!")
        for n in range(min(5, nd)):
            print(f"    lag {n}: MATLAB={c_mat0[n]:.8e}  BMS={c_bms0[n]:.8e}")

    # ---- Motion-blur-only comparison (sigma=0) --------------------------------
    print()
    print("─" * 70)
    print("  Noise-free (sigma=0): isolate motion-blur term")
    c_mat_nb, c_bms_nb = first_row_comparison(K, alpha, 0.0, dt, frac, nd)
    reldiff_nb = (c_bms_nb - c_mat_nb) / np.maximum(np.abs(c_mat_nb), eps) * 100
    print(f"  {'lag (n)':<8}  {'MATLAB c(n)':>14}  {'bayesmsd c(n)':>14}  {'rel diff %':>10}")
    print("  " + "-" * 52)
    for n in range(nd):
        print(f"  {n:<8}  {c_mat_nb[n]:>14.6e}  {c_bms_nb[n]:>14.6e}  {reldiff_nb[n]:>10.4f}")

    # ---- Cross-fit (optional): fit bayesmsd model to MATLAB covariance --------
    if args.cross_fit:
        print()
        print("─" * 70)
        print("  Cross-fit: fitting bayesmsd model to MATLAB covariance target")
        print("  (shows alpha shift purely from covariance model difference)")
        K_fit, a_fit, s_fit, resid = fit_alpha_from_cov(
            c_mat, K, alpha, sigma, dt, frac, nd, method='bms')
        print(f"  MATLAB target alpha : {alpha:.4f}  K={K:.5f}  sigma={sigma:.5f}")
        print(f"  bayesmsd fit alpha  : {a_fit:.4f}  K={K_fit:.5f}  sigma={s_fit:.5f}")
        print(f"  Delta alpha         : {a_fit - alpha:+.4f}")
        print(f"  Residual norm       : {resid:.4e}")

        print()
        print("  Cross-fit: fitting MATLAB model to bayesmsd covariance target")
        K_fit2, a_fit2, s_fit2, resid2 = fit_alpha_from_cov(
            c_bms, K, alpha, sigma, dt, frac, nd, method='matlab')
        print(f"  bayesmsd target alpha : {alpha:.4f}  K={K:.5f}  sigma={sigma:.5f}")
        print(f"  MATLAB fit alpha      : {a_fit2:.4f}  K={K_fit2:.5f}  sigma={s_fit2:.5f}")
        print(f"  Delta alpha           : {a_fit2 - alpha:+.4f}")
        print(f"  Residual norm         : {resid2:.4e}")

    print()
    print("=" * 70)
    print("  Summary")
    print("─" * 70)
    print(f"  With motion blur f={frac}, alpha={alpha}, nd={nd}:")
    print(f"  Max relative difference in covariance: {max_rel:.4f}%")
    if max_rel < 0.01:
        print("  → The two covariance models are NUMERICALLY IDENTICAL.")
        print("  → The fitted alpha discrepancy is NOT due to a kernel difference.")
        print("  → Likely causes: optimizer landscape, parametrization, or")
        print("    different sufficient-statistic handling (per-track vs pooled).")
    else:
        print("  → The two covariance models DIFFER at the above level.")
        print("  → This kernel mismatch can systematically shift fitted alpha.")
    print()


if __name__ == "__main__":
    main()

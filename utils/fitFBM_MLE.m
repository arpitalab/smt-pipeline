function [params, loglik, CI] = fitFBM_MLE(tracks, dt, varargin)
%FITFBM_MLE  Maximum-likelihood estimation of fBM parameters from 2D tracks.
%
%   [params, loglik]     = fitFBM_MLE(tracks, dt)
%   [params, loglik, CI] = fitFBM_MLE(tracks, dt, 'CIMethod', 'profile', ...)
%   [params, loglik, CI] = fitFBM_MLE(tracks, dt, 'CIMethod', 'mcmc', ...)
%
%   Fits the generalized diffusion coefficient K (µm²/s^alpha), anomalous
%   exponent alpha, and localization precision sigma (µm) by maximising the
%   multivariate-Gaussian log-likelihood of the 2D displacement time series.
%
%   COVARIANCE MODEL (Backlund et al. Phys. Rev. E 2015)
%   -------------------------------------------------------
%   For a subtrack of length L, the (L-1) consecutive displacements are
%   jointly Gaussian with Toeplitz covariance.  The n-th off-diagonal is:
%
%     Sigma(n) = K * [Psi((n+1)*dt,Te,a) + Psi(|n-1|*dt,Te,a) - 2*Psi(n*dt,Te,a)]
%                + 2*sigma^2*(n==0)  -  sigma^2*(n==1)
%
%   where the exposure-time-averaged fBM kernel is:
%
%     Psi(tau,Te,a) = [(tau+Te)^(a+2) + |tau-Te|^(a+2) - 2*tau^(a+2)]
%                    / [Te^2 * (a+1) * (a+2)]
%
%   and Te = ExposureFraction * dt.  For stroboscopic imaging (Te→0),
%   Psi(tau,Te,a) → tau^a, recovering the standard fBM Toeplitz structure.
%   x and y are treated as independent (isotropic diffusion).
%
%   PARAMETRIZATION
%   -------------------------------------------------------
%   Internally, optimisation runs in unconstrained space:
%     theta(1) = log(K * alpha)    the "transport parameter" Ka (µm²/s^alpha)
%     theta(2) = logit(alpha/2)    exponent mapped to (-inf, inf)
%     theta(3) = log(sigma)        localization precision
%
%   Fitting log(K*alpha) rather than log(K) directly reduces the posterior
%   correlation between K and alpha (K*alpha appears naturally in the Fisher
%   information cross-term).  The back-transformation is:
%     alpha = 2 * sigmoid(theta(2))
%     K     = exp(theta(1)) / alpha
%     sigma = exp(theta(3))
%
%   CONFIDENCE INTERVALS
%   -------------------------------------------------------
%   'profile'  (default) — profile log-likelihood for each parameter,
%              maximising over the other two.  Gives likelihood-ratio CI
%              that correctly handles asymmetry and banana-shaped posteriors.
%              95% threshold: MLE loglik - chi2_95/2 = MLE - 1.9208.
%
%   'mcmc'     — adaptive Metropolis-Hastings sampling of the posterior
%              (flat prior → posterior ∝ likelihood).  Returns full chains
%              in CI.samples plus 2.5/97.5 percentile bounds.
%
%   'hessian'  — Laplace approximation via numerical Hessian.  Fast but
%              poor for this problem; use only for quick diagnostics.
%
%   Inputs:
%     tracks  - Cell array of trajectory matrices, each Nx2 or larger
%               (col 1 = x [µm], col 2 = y [µm], contiguous frames).
%     dt      - Frame interval in seconds.
%
%   Name-value options:
%     'ExposureFraction'  - Te/dt ratio (0 = stroboscopic, default: 1).
%     'SubtrackLength'    - Target frames per subtrack (default: 20).
%     'MinSubtrackLength' - Minimum frames to include a subtrack (default: 10).
%     'InitialGuess'      - [K0, alpha0, sigma0] starting point.
%     'CIMethod'          - 'profile' | 'profile+mcmc' | 'mcmc' | 'hessian' | 'none'
%                           (default: 'profile').
%                           'profile+mcmc' runs profiles first, uses the profile-derived
%                           marginal widths to calibrate the MCMC proposal, and defaults
%                           to a short burn-in (300 steps) since the chain starts at the
%                           MLE and the proposal is already well-scaled.
%     'ProfilePoints'     - Grid points per parameter for profile CI (default: 60).
%     'MCMCSamples'       - Posterior samples after burn-in (default: 5000).
%     'MCMCBurnin'        - MCMC burn-in steps (default: 2000).
%     'Verbose'           - Print optimisation progress (default: false).
%
%   Outputs:
%     params  - Struct:
%                 .K              generalized diffusion coeff (µm²/s^alpha)
%                 .Ka             K * alpha  (the fitted transport parameter)
%                 .alpha          anomalous exponent in (0, 2)
%                 .sigma          localization precision (µm)
%                 .loglik         total log-likelihood at MLE
%                 .n_subtracks    number of subtracks used
%                 .n_displacements  total scalar displacement observations
%     loglik  - Scalar log-likelihood at MLE (= params.loglik).
%     CI      - Struct with 95% bounds:
%                 .K, .Ka, .alpha, .sigma  each [lower, upper]
%                 .method                  string identifying CI method
%                 .samples   (mcmc only)   Nx3 matrix [K, alpha, sigma]
%               Empty struct if CIMethod = 'none'.

p = inputParser;
addRequired(p,  'tracks',                  @iscell);
addRequired(p,  'dt',                      @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ExposureFraction',   1,   @isnumeric);
addParameter(p, 'SubtrackLength',    20,   @isnumeric);
addParameter(p, 'MinSubtrackLength', 10,   @isnumeric);
addParameter(p, 'InitialGuess',      [],   @isnumeric);
addParameter(p, 'CIMethod',    'profile',  @(x) ismember(lower(x), {'profile','profile+mcmc','mcmc','hessian','none'}));
addParameter(p, 'ProfilePoints',     60,   @isnumeric);
addParameter(p, 'MCMCSamples',     5000,   @isnumeric);
addParameter(p, 'MCMCBurnin',      2000,   @isnumeric);
addParameter(p, 'Verbose',        false,   @islogical);
parse(p, tracks, dt, varargin{:});
o = p.Results;

frac = o.ExposureFraction;
Te   = frac * dt;
Lsub = round(o.SubtrackLength);
Lmin = round(o.MinSubtrackLength);

if Lmin < 3
    error('fitFBM_MLE: MinSubtrackLength must be at least 3.');
end
if Lsub < Lmin
    error('fitFBM_MLE: SubtrackLength must be >= MinSubtrackLength.');
end

% -------------------------------------------------------------------------
% 1.  Collect subtracks and precompute sufficient statistics
% -------------------------------------------------------------------------
% The Gaussian log-likelihood depends on the raw displacements only through
% the n×n sample covariance matrix  S_L = DX_L'*DX_L + DY_L'*DY_L
% (where DX_L is N_L×n_L, n_L = L-1).  Computing S_L once makes each
% subsequent likelihood evaluation O(n_L^3) regardless of N_L — a
% ~N_L/n_L fold speedup over looping over subtracks, which matters when
% N_L is in the tens of thousands.

Smat   = cell(Lsub, 1);   % Smat{L}  : (L-1)×(L-1) sufficient-statistic matrix
N_subs = zeros(Lsub, 1);  % N_subs(L): number of subtracks of length L

for k = 1:numel(tracks)
    tr = tracks{k};
    if isempty(tr) || size(tr, 1) < Lmin
        continue;
    end
    x = tr(:, 1);
    y = tr(:, 2);
    n = size(tr, 1);
    i0 = 1;
    while i0 <= n
        i1 = min(i0 + Lsub - 1, n);
        L  = i1 - i0 + 1;
        if L >= Lmin
            dx = diff(x(i0:i1));   % (L-1)×1
            dy = diff(y(i0:i1));
            if isempty(Smat{L})
                Smat{L} = dx*dx' + dy*dy';
            else
                Smat{L} = Smat{L} + dx*dx' + dy*dy';   % accumulate outer products
            end
            N_subs(L) = N_subs(L) + 1;
        end
        i0 = i0 + Lsub;
    end
end

nSubtracks = sum(N_subs);
nDisp      = 2 * sum(N_subs .* max(0, (1:Lsub)' - 1));

if nSubtracks == 0
    error('fitFBM_MLE: no subtracks of length >= %d found.', Lmin);
end

% -------------------------------------------------------------------------
% 2.  Maximum likelihood optimisation (unconstrained, theta-space)
% -------------------------------------------------------------------------
if isempty(o.InitialGuess)
    K0 = 1e-3; alpha0 = 0.8; sigma0 = 0.05;
else
    K0 = o.InitialGuess(1); alpha0 = o.InitialGuess(2); sigma0 = o.InitialGuess(3);
end

theta0 = encode_theta(K0, alpha0, sigma0);

negll = @(th) -compute_loglik(th, Smat, N_subs, Lsub, dt, Te);

opts_fmin = optimoptions('fminunc', ...
    'MaxFunctionEvaluations', 8000, ...
    'MaxIterations',          3000, ...
    'OptimalityTolerance',    1e-9, ...
    'StepTolerance',          1e-11, ...
    'Display', local_ifelse(o.Verbose, 'iter', 'none'));

[theta_hat, nll_hat] = fminunc(negll, theta0, opts_fmin);

[K_hat, alpha_hat, sigma_hat] = decode_theta(theta_hat);
Ka_hat = K_hat * alpha_hat;   % = exp(theta_hat(1))

loglik = -nll_hat;

params.K               = K_hat;
params.Ka              = Ka_hat;
params.alpha           = alpha_hat;
params.sigma           = sigma_hat;
params.loglik          = loglik;
params.n_subtracks     = nSubtracks;
params.n_displacements = nDisp;

CI = struct('method', lower(o.CIMethod));

if strcmp(lower(o.CIMethod), 'none')
    return;
end

ci_method  = lower(o.CIMethod);
do_profile = ismember(ci_method, {'profile', 'profile+mcmc'});
do_mcmc    = ismember(ci_method, {'mcmc',    'profile+mcmc'});
do_hessian = strcmp(ci_method, 'hessian');

% Always compute Hessian — used for proposal covariance and search-range estimates.
H     = numerical_hessian(negll, theta_hat);
invH  = safe_inv(H);
se_th = sqrt(max(1e-10, diag(invH)));   % Hessian-based marginal SEs in theta-space

% -------------------------------------------------------------------------
% 3a.  Profile log-likelihood (computed for 'profile' and 'profile+mcmc')
% -------------------------------------------------------------------------
%
%  Strategy: for each parameter theta_i, sweep a grid while maximising the
%  log-likelihood over the other two (conditional optimisation).  Grid points
%  are swept outward from the MLE so that each optimisation warm-starts from
%  the previous solution — this is critical for convergence near the edges.
%  Grid range: MLE ± 5 Hessian-SEs, wide enough to bracket the 95% CI.
%
%  The 95% likelihood-ratio threshold is  MLE_ll - chi2_95/2 = MLE - 1.9208.
%  CI bounds are where the profile crosses this threshold.

if do_profile
    chi2_thresh = loglik - 1.9208;

    opts_prof = optimoptions('fminunc', ...
        'MaxFunctionEvaluations', 3000, ...
        'OptimalityTolerance',    1e-8, ...
        'Display', 'none');

    n_grid    = o.ProfilePoints;
    th_ranges = [theta_hat' - 5*se_th, theta_hat' + 5*se_th];

    for i_param = 1:3
        th_grid   = linspace(th_ranges(i_param,1), th_ranges(i_param,2), n_grid);
        pll       = nan(1, n_grid);
        idx_other = setdiff(1:3, i_param);

        [~, i_mle] = min(abs(th_grid - theta_hat(i_param)));

        % Sweep upward from MLE, warm-starting each step from the previous optimum
        th_warm = theta_hat(idx_other);
        for i = i_mle : n_grid
            [pll(i), th_warm] = profile_step(th_grid(i), i_param, th_warm, negll, opts_prof);
        end
        % Sweep downward from MLE, resetting warm-start at MLE
        th_warm = theta_hat(idx_other);
        for i = i_mle-1 : -1 : 1
            [pll(i), th_warm] = profile_step(th_grid(i), i_param, th_warm, negll, opts_prof);
        end

        above  = pll >= chi2_thresh;
        lo_idx = find(above, 1, 'first');
        hi_idx = find(above, 1, 'last');

        CI.profiles(i_param).theta_grid = th_grid;
        CI.profiles(i_param).pll        = pll;
        CI.theta_bounds(i_param, :)     = [th_grid(lo_idx), th_grid(hi_idx)];
    end

    CI.Ka    = sort([exp(CI.theta_bounds(1,1)),  exp(CI.theta_bounds(1,2))]);
    CI.alpha = sort([decode_alpha(CI.theta_bounds(2,1)), decode_alpha(CI.theta_bounds(2,2))]);
    CI.sigma = sort([exp(CI.theta_bounds(3,1)),  exp(CI.theta_bounds(3,2))]);
    % K CI: K = Ka/alpha, so extremes come from combining Ka and alpha profile bounds
    CI.K     = sort([CI.Ka(1)/CI.alpha(2), CI.Ka(2)/CI.alpha(1)]);
end

% -------------------------------------------------------------------------
% 3b.  Adaptive Metropolis-Hastings MCMC
% -------------------------------------------------------------------------
%
%  Starting point: always the MLE (avoids the expensive "find the mode"
%  phase of burn-in; decorrelation burn-in is all that remains).
%
%  Proposal covariance:
%    'mcmc' alone  — Hessian-based covariance (standard Roberts & Rosenthal
%                    (2.38²/d) scaling).
%    'profile+mcmc' — Profile-calibrated covariance: the Hessian correlation
%                    structure is preserved but each marginal is rescaled to
%                    the profile-derived standard error, which is accurate
%                    even for non-Gaussian (e.g. skewed) posteriors.
%                    Burn-in defaults to 300 steps instead of 2000.
%
%  Adaptation: proposal covariance is updated every 300 steps during burn-in
%  from the empirical chain covariance, so it converges to the true posterior
%  covariance regardless of the initial proposal.

if do_mcmc
    d         = 3;
    n_samples = o.MCMCSamples;
    n_burnin  = o.MCMCBurnin;

    % For profile+mcmc the default burn-in is much shorter (user can override
    % via 'MCMCBurnin').  Only apply shortcut when MCMCBurnin was not
    % explicitly set (i.e. it is still at the parser default of 2000).
    if strcmp(ci_method, 'profile+mcmc') && n_burnin == 2000
        n_burnin = 300;
    end

    if do_profile
        % Construct proposal covariance: profile-accurate marginals + Hessian correlations.
        % profile_se = 95% CI half-width / 1.96  (marginal SE in theta-space)
        profile_se = (CI.theta_bounds(:,2) - CI.theta_bounds(:,1)) / (2 * 1.96);
        profile_se = max(profile_se, se_th * 0.1);   % guard against collapsed CI

        % Correlation matrix from Hessian
        d_hess  = sqrt(diag(invH));
        corr_H  = invH ./ (d_hess * d_hess');
        corr_H  = (corr_H + corr_H') / 2;            % enforce symmetry

        % Scale Hessian correlations by profile-derived marginals
        prop_cov = (2.38^2 / d) * diag(profile_se) * corr_H * diag(profile_se);
    else
        prop_cov = (2.38^2 / d) * invH;
    end

    [L_prop, flag] = chol(prop_cov + 1e-8*eye(d), 'lower');
    if flag ~= 0
        L_prop = diag(se_th * 2.38 / sqrt(d));   % diagonal fallback
    end

    n_total    = n_burnin + n_samples;
    chain      = zeros(n_total, d);
    chain(1,:) = theta_hat;
    ll_cur     = loglik;
    n_accept   = 0;

    for i = 2:n_total
        theta_prop = chain(i-1,:)' + L_prop * randn(d, 1);
        ll_prop    = -negll(theta_prop');

        if log(rand) < ll_prop - ll_cur
            chain(i,:) = theta_prop';
            ll_cur     = ll_prop;
            n_accept   = n_accept + 1;
        else
            chain(i,:) = chain(i-1,:);
        end

        % Empirical adaptation during burn-in (every 300 steps, using last 600)
        if mod(i, 300) == 0 && i <= n_burnin
            recent = chain(max(1, i-600):i, :);
            if size(recent, 1) > d + 1
                C = cov(recent) * (2.38^2 / d);
                [L_new, ok] = chol(C + 1e-8*eye(d), 'lower');
                if ok == 0
                    L_prop = L_new;
                end
            end
        end
    end

    accept_rate = n_accept / n_total;
    if o.Verbose
        fprintf('MCMC: %d samples, %d burn-in, acceptance %.1f%%\n', ...
            n_samples, n_burnin, 100 * accept_rate);
    end
    if accept_rate < 0.10 || accept_rate > 0.70
        warning('fitFBM_MLE:mcmcAcceptance', ...
            'MCMC acceptance rate %.1f%% is outside the recommended 10-60%% range.', ...
            100 * accept_rate);
    end

    post = chain(n_burnin+1:end, :);

    K_samp     = zeros(n_samples, 1);
    alpha_samp = zeros(n_samples, 1);
    sigma_samp = zeros(n_samples, 1);
    Ka_samp    = zeros(n_samples, 1);
    for i = 1:n_samples
        [K_samp(i), alpha_samp(i), sigma_samp(i)] = decode_theta(post(i,:));
        Ka_samp(i) = K_samp(i) * alpha_samp(i);
    end

    CI.K           = quantile(K_samp,     [0.025, 0.975]);
    CI.Ka          = quantile(Ka_samp,    [0.025, 0.975]);
    CI.alpha       = quantile(alpha_samp, [0.025, 0.975]);
    CI.sigma       = quantile(sigma_samp, [0.025, 0.975]);
    CI.samples     = [K_samp, alpha_samp, sigma_samp, Ka_samp];
    CI.accept_rate = accept_rate;
    CI.n_burnin    = n_burnin;
end

% -------------------------------------------------------------------------
% 3c.  Hessian / Laplace approximation (fast diagnostic only)
% -------------------------------------------------------------------------
if do_hessian
    se = sqrt(max(0, diag(invH)));
    z  = 1.96;
    lo = theta_hat - z * se';
    hi = theta_hat + z * se';

    CI.Ka    = sort([exp(lo(1)), exp(hi(1))]);
    CI.alpha = sort([decode_alpha(lo(2)), decode_alpha(hi(2))]);
    CI.sigma = sort([exp(lo(3)), exp(hi(3))]);
    CI.K     = sort([CI.Ka(1)/CI.alpha(2), CI.Ka(2)/CI.alpha(1)]);
end

end % fitFBM_MLE


% =========================================================================
%  Local functions
% =========================================================================

function ll = compute_loglik(theta, Smat, N_subs, Lsub, dt, Te)
%COMPUTE_LOGLIK  Total log-likelihood via precomputed sufficient statistics.
%
%   The Gaussian likelihood for N_L subtracks of length L with displacements
%   pooled from x and y depends on the raw data only through the n×n matrix
%     S_L = sum_k (dx_k*dx_k' + dy_k*dy_k')   [n = L-1]
%   so each evaluation costs O(n^3) regardless of how many subtracks there are.
%
%   ℓ = Σ_L  N_L*(-n*log(2π) - log|Σ_L|) - (1/2)*tr(Σ_L^{-1} * S_L)

[K, alpha, sigma] = decode_theta(theta);
if K <= 0 || alpha <= 0 || alpha >= 2 || sigma < 0 || ~isreal([K, alpha, sigma])
    ll = -1e15;
    return;
end

ll = 0;
for L = 2:Lsub
    if N_subs(L) == 0
        continue;
    end
    nd   = L - 1;
    NL   = N_subs(L);

    Sigma = build_cov(nd, K, alpha, sigma, dt, Te);

    [R, flag] = chol(Sigma);
    if flag ~= 0
        ll = -1e15;
        return;
    end

    log_det = 2 * sum(log(diag(R)));

    % tr(Sigma^{-1} * S) = tr(R^{-1} * (R'^{-1} * S))
    % Computed via two triangular solves on the n×n matrix S_L
    V  = R' \ Smat{L};   % R'*V = S  →  V = (R')^{-1} * S
    qf = trace(R \ V);   % R*W = V   →  trace(W) = tr(Sigma^{-1}*S)

    ll = ll + NL * (-nd * log(2*pi) - log_det) - 0.5 * qf;
end
end


function Sigma = build_cov(nd, K, alpha, sigma, dt, Te)
%BUILD_COV  Build the nd×nd Toeplitz displacement covariance matrix.

m = (0 : nd-1)';
psi_p = arrayfun(@(m_) psi_fbm((m_+1)*dt, Te, alpha), m);
psi_0 = arrayfun(@(m_) psi_fbm( m_   *dt, Te, alpha), m);
psi_m = arrayfun(@(m_) psi_fbm(abs(m_-1)*dt, Te, alpha), m);

c    = K * (psi_p + psi_m - 2*psi_0);
c(1) = c(1) + 2*sigma^2;      % diagonal (m = 0)
if nd > 1
    c(2) = c(2) - sigma^2;    % nearest-neighbor anti-correlation (m = 1)
end

Sigma = toeplitz(c);
end


function val = psi_fbm(tau, Te, alpha)
%PSI_FBM  Exposure-time-averaged fBM covariance kernel.
%
%   Derived from the double time-average of the fBM covariance over the
%   camera exposure window [t, t+Te]:
%
%     Psi(tau, Te, alpha) = [(tau+Te)^(a+2) + |tau-Te|^(a+2) - 2*tau^(a+2)]
%                           / [Te^2 * (a+1) * (a+2)]
%
%   Valid for all tau >= 0.  Limit Te -> 0 gives tau^alpha (stroboscopic).

if Te < 1e-10
    val = tau^alpha;
else
    val = ((tau + Te)^(alpha+2) + abs(tau - Te)^(alpha+2) - 2*tau^(alpha+2)) ...
          / (Te^2 * (alpha+1) * (alpha+2));
end
end


function theta = encode_theta(K, alpha, sigma)
%ENCODE_THETA  Map (K, alpha, sigma) → unconstrained theta.
%   theta(1) = log(K * alpha)    transport parameter Ka
%   theta(2) = logit(alpha/2)    exponent on (-inf, inf)
%   theta(3) = log(sigma)
theta = [log(K * alpha), log((alpha/2) / (1 - alpha/2)), log(sigma)];
end


function [K, alpha, sigma] = decode_theta(theta)
%DECODE_THETA  Map unconstrained theta → (K, alpha, sigma).
alpha = 2 / (1 + exp(-theta(2)));     % sigmoid: alpha in (0, 2)
K     = exp(theta(1)) / alpha;        % K = Ka / alpha
sigma = exp(theta(3));
end


function alpha = decode_alpha(th2)
%DECODE_ALPHA  Scalar helper for alpha from theta(2).
alpha = 2 / (1 + exp(-th2));
end


function [pll_i, th_other_opt] = profile_step(th_fixed, i_param, th_other0, negll, opts)
%PROFILE_STEP  One grid-point evaluation for profile likelihood.
%   Returns th_other_opt for warm-starting the next grid point.
idx_other = setdiff(1:3, i_param);
neg_ll_cond = @(th_o) negll(insert_param(th_o, th_fixed, i_param, idx_other));
[th_other_opt, nll] = fminunc(neg_ll_cond, th_other0, opts);
pll_i = -nll;
end


function theta = insert_param(th_other, th_fixed, i_fix, idx_other)
%INSERT_PARAM  Reconstruct full theta from conditional and fixed parts.
theta = zeros(1, 3);
theta(i_fix)    = th_fixed;
theta(idx_other) = th_other;
end


function H = numerical_hessian(f, x0)
%NUMERICAL_HESSIAN  Central-difference Hessian of scalar f at x0.
n  = numel(x0);
H  = zeros(n);
ep = 1e-4;
f0 = f(x0);
for i = 1:n
    for j = i:n
        ei = zeros(1, n); ei(i) = ep;
        ej = zeros(1, n); ej(j) = ep;
        if i == j
            H(i,i) = (f(x0+ei) - 2*f0 + f(x0-ei)) / ep^2;
        else
            H(i,j) = (f(x0+ei+ej) - f(x0+ei-ej) - f(x0-ei+ej) + f(x0-ei-ej)) / (4*ep^2);
            H(j,i) = H(i,j);
        end
    end
end
end


function invH = safe_inv(H)
%SAFE_INV  Regularised inversion: add small diagonal if H is not PD.
[~, flag] = chol(H);
if flag ~= 0
    lam_min = min(eig(H));
    H = H + (abs(lam_min) + 1e-6) * eye(size(H, 1));
end
invH = inv(H);
end


function s = local_ifelse(cond, t, f)
if cond; s = t; else; s = f; end
end

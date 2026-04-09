function results = fitFBM_pertracks(tracks, dt, varargin)
%FITFBM_PERTRACKS  Per-track fBM parameter estimation (K and alpha distribution).
%
%   results = fitFBM_pertracks(tracks, dt)
%   results = fitFBM_pertracks(tracks, dt, 'SigmaFixed', 0.03, ...)
%
%   Fits (K, alpha) independently to each track by maximising the
%   multivariate-Gaussian log-likelihood of the displacement time series.
%   Sigma (localization precision) is either fixed to a known value or fitted
%   jointly (3-parameter problem per track).
%
%   Fixing sigma is strongly recommended when tracks are short (< 30 frames):
%   with only ~2*(L-1) scalar observations and 3 parameters, jointly-fitted
%   sigma is poorly identified per track.  Obtain sigma once from the
%   population MLE (fitFBM_MLE / getFBMParameters) and fix it here.
%
%   Inputs:
%     tracks    - Cell array of trajectory matrices (col 1=x µm, col 2=y µm).
%     dt        - Frame interval (s).
%
%   Name-value options:
%     'ExposureFraction'   - Te/dt (default: 1).
%     'SigmaFixed'         - Scalar sigma (µm) to hold fixed.  NaN = fit jointly.
%     'MinTrackLength'     - Minimum frames per track (default: 15).
%     'MaxSubtrackLength'  - Maximum frames used per fit.  Tracks longer than
%                            this are split into non-overlapping subtracks of
%                            this length; their sufficient statistics S are
%                            accumulated before fitting, so all data contribute
%                            but the Toeplitz matrix stays MaxSubtrackLength-1
%                            square.  This keeps the lag range consistent with
%                            fitFBM_MLE (set to the same SubtrackLength, e.g. 20).
%                            Default: 0 = use the full track as one subtrack.
%     'MinStepVar'         - Minimum mean squared step size (µm²) to include a
%                            track.  Tracks below this threshold are essentially
%                            immobile; they contribute near-zero displacements
%                            whose Toeplitz fit collapses alpha toward zero.
%                            Set to e.g. 4*sigma^2 to exclude stuck particles.
%                            Default: 0 (no filter).
%     'InitK'              - Initial K guess (µm²/s^a, default: 1e-3).
%     'InitAlpha'          - Initial alpha guess (default: 0.8).
%     'InitSigma'          - Initial sigma guess when fitting jointly (default: 0.05).
%     'Verbose'            - Print progress every 1000 tracks (default: false).
%
%   Output struct fields:
%     .alpha        Nx1  per-track anomalous exponent
%     .K            Nx1  per-track generalized diffusion coeff (µm²/s^alpha)
%     .Ka           Nx1  K * alpha  (the fitted transport parameter)
%     .sigma        Nx1  per-track sigma (= SigmaFixed for all if fixed)
%     .track_length Nx1  number of frames in each track
%     .loglik       Nx1  per-track log-likelihood at optimum
%     .converged    Nx1  logical: true if fminunc exited cleanly
%     .sigma_fixed  scalar fixed sigma (NaN if fitted jointly)
%     .dt           frame interval used
%     .frac         ExposureFraction used

p = inputParser;
addRequired(p,  'tracks',                @iscell);
addRequired(p,  'dt',                    @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ExposureFraction',    1,  @isnumeric);
addParameter(p, 'SigmaFixed',        NaN,  @isnumeric);
addParameter(p, 'MinTrackLength',     15,  @isnumeric);
addParameter(p, 'MaxSubtrackLength',   0,  @isnumeric);
addParameter(p, 'MinStepVar',          0,  @isnumeric);
addParameter(p, 'InitK',           1e-3,  @isnumeric);
addParameter(p, 'InitAlpha',        0.8,  @isnumeric);
addParameter(p, 'InitSigma',       0.05,  @isnumeric);
addParameter(p, 'Verbose',         false, @islogical);
parse(p, tracks, dt, varargin{:});
o = p.Results;

frac    = o.ExposureFraction;
Te      = frac * dt;
Lmin    = round(o.MinTrackLength);
Lsub    = round(o.MaxSubtrackLength);   % 0 means use full track
fix_sig = ~isnan(o.SigmaFixed);
sigma0  = local_ifelse(fix_sig, o.SigmaFixed, o.InitSigma);

% Shared fminunc options — created once, reused for every track
opts = optimoptions('fminunc', ...
    'MaxFunctionEvaluations', 600, ...
    'OptimalityTolerance',    1e-7, ...
    'StepTolerance',          1e-9, ...
    'Display',                'none');

% Pre-filter tracks by length
keep = cellfun(@(t) size(t,1) >= Lmin, tracks);
tracks = tracks(keep);

% Pre-filter by minimum step variance (removes essentially immobile particles)
if o.MinStepVar > 0
    keep2 = cellfun(@(t) mean(diff(t(:,1)).^2 + diff(t(:,2)).^2) >= o.MinStepVar, tracks);
    tracks = tracks(keep2);
end

N = numel(tracks);

if N == 0
    error('fitFBM_pertracks: no tracks >= MinTrackLength=%d frames.', Lmin);
end

% Pre-allocate outputs
alpha_out  = nan(N, 1);
K_out      = nan(N, 1);
Ka_out     = nan(N, 1);
sigma_out  = nan(N, 1);
ll_out     = nan(N, 1);
len_out    = zeros(N, 1);
conv_out   = false(N, 1);

% Initial theta (shared warm-start; each track updates its own)
if fix_sig
    % 2D: theta = [log(K*alpha), logit(alpha/2)]
    th0_2d = [log(o.InitK * o.InitAlpha), log((o.InitAlpha/2)/(1 - o.InitAlpha/2))];
else
    % 3D: theta = [log(K*alpha), logit(alpha/2), log(sigma)]
    th0_3d = [log(o.InitK * o.InitAlpha), log((o.InitAlpha/2)/(1 - o.InitAlpha/2)), log(sigma0)];
end

for k = 1:N
    tr = tracks{k};
    x  = tr(:, 1);
    y  = tr(:, 2);
    L  = size(tr, 1);

    % Build sufficient statistic S = Σ(Δx·Δxᵀ + Δy·Δyᵀ).
    % If MaxSubtrackLength is set and the track is longer, split into
    % non-overlapping subtracks of that length and accumulate their S
    % matrices.  This caps the Toeplitz matrix size (= lag range) to be
    % consistent with fitFBM_MLE while using all available frames.
    if Lsub > 0 && L > Lsub
        nd    = Lsub - 1;
        n_sub = floor(L / Lsub);
        S     = zeros(nd, nd);
        for isub = 1:n_sub
            i1  = (isub-1)*Lsub + 1;
            i2  = i1 + Lsub - 1;
            dxi = diff(x(i1:i2));
            dyi = diff(y(i1:i2));
            S   = S + dxi*dxi' + dyi*dyi';
        end
    else
        nd    = L - 1;
        n_sub = 1;
        dx    = diff(x);
        dy    = diff(y);
        S     = dx*dx' + dy*dy';
    end

    if fix_sig
        neg_ll = @(th) -track_loglik_2d(th, S, nd, n_sub, dt, Te, o.SigmaFixed);
        [th_hat, nll, ~, out] = fminunc(neg_ll, th0_2d, opts);
        alpha = 2 / (1 + exp(-th_hat(2)));
        K     = exp(th_hat(1)) / alpha;
        sig   = o.SigmaFixed;
    else
        neg_ll = @(th) -track_loglik_3d(th, S, nd, n_sub, dt, Te);
        [th_hat, nll, ~, out] = fminunc(neg_ll, th0_3d, opts);
        alpha = 2 / (1 + exp(-th_hat(2)));
        K     = exp(th_hat(1)) / alpha;
        sig   = exp(th_hat(3));
    end

    alpha_out(k) = alpha;
    K_out(k)     = K;
    Ka_out(k)    = K * alpha;
    sigma_out(k) = sig;
    ll_out(k)    = -nll;
    len_out(k)   = L;
    conv_out(k)  = (out.firstorderopt < 1e-4);

    if o.Verbose && mod(k, 1000) == 0
        fprintf('  fitFBM_pertracks: %d / %d tracks\n', k, N);
    end
end

results.alpha        = alpha_out;
results.K            = K_out;
results.Ka           = Ka_out;
results.sigma        = sigma_out;
results.track_length = len_out;
results.loglik       = ll_out;
results.converged    = conv_out;
results.sigma_fixed  = local_ifelse(fix_sig, o.SigmaFixed, NaN);
results.dt           = dt;
results.frac         = frac;
results.n_tracks     = N;

if o.Verbose
    fprintf('fitFBM_pertracks: done. median alpha=%.3f, median K=%.4g\n', ...
        median(alpha_out, 'omitnan'), median(K_out, 'omitnan'));
end
end


% =========================================================================
%  Local functions (self-contained; replicate the covariance model from
%  fitFBM_MLE so this file has no dependency on it)
% =========================================================================

function ll = track_loglik_2d(theta2, S, nd, n_sub, dt, Te, sigma)
%TRACK_LOGLIK_2D  Log-likelihood with fixed sigma; theta = [log(Ka), logit(a/2)].
alpha = 2 / (1 + exp(-theta2(2)));
K     = exp(theta2(1)) / alpha;
if K <= 0 || alpha <= 0 || alpha >= 2
    ll = -1e15; return;
end
ll = gauss_loglik(S, nd, n_sub, K, alpha, sigma, dt, Te);
end


function ll = track_loglik_3d(theta3, S, nd, n_sub, dt, Te)
%TRACK_LOGLIK_3D  Log-likelihood; theta = [log(Ka), logit(a/2), log(sigma)].
alpha = 2 / (1 + exp(-theta3(2)));
K     = exp(theta3(1)) / alpha;
sigma = exp(theta3(3));
if K <= 0 || alpha <= 0 || alpha >= 2 || sigma < 0
    ll = -1e15; return;
end
ll = gauss_loglik(S, nd, n_sub, K, alpha, sigma, dt, Te);
end


function ll = gauss_loglik(S, nd, n_sub, K, alpha, sigma, dt, Te)
%GAUSS_LOGLIK  Pooled log-likelihood via sufficient statistic S.
%
%   For n_sub non-overlapping subtracks of length nd, each contributing an
%   (nd x nd) outer-product matrix, S = sum of all subtrack outer products.
%   The combined log-likelihood is:
%     ll = n_sub * (-nd*log(2pi) - log|Sigma|) - (1/2)*tr(Sigma^{-1}*S)
%
%   For a single subtrack (n_sub=1) this reduces to the standard formula.
Sigma = build_cov(nd, K, alpha, sigma, dt, Te);
[R, flag] = chol(Sigma);
if flag ~= 0
    ll = -1e15; return;
end
log_det = 2 * sum(log(diag(R)));
V  = R' \ S;
qf = trace(R \ V);
ll = n_sub * (-nd * log(2*pi) - log_det) - 0.5 * qf;
end


function Sigma = build_cov(nd, K, alpha, sigma, dt, Te)
%BUILD_COV  Toeplitz displacement covariance matrix.
m    = (0:nd-1)';
pp   = arrayfun(@(m_) psi_fbm((m_+1)*dt, Te, alpha), m);
p0   = arrayfun(@(m_) psi_fbm( m_   *dt, Te, alpha), m);
pm   = arrayfun(@(m_) psi_fbm(abs(m_-1)*dt, Te, alpha), m);
c    = K * (pp + pm - 2*p0);
c(1) = c(1) + 2*sigma^2;
if nd > 1
    c(2) = c(2) - sigma^2;
end
Sigma = toeplitz(c);
end


function val = psi_fbm(tau, Te, alpha)
%PSI_FBM  Exposure-time-averaged fBM covariance kernel.
if Te < 1e-10
    val = tau^alpha;
else
    val = ((tau+Te)^(alpha+2) + abs(tau-Te)^(alpha+2) - 2*tau^(alpha+2)) ...
          / (Te^2 * (alpha+1) * (alpha+2));
end
end


function s = local_ifelse(cond, t, f)
if cond; s = t; else; s = f; end
end

function out = compareStateArrayPosteriors(boot_A, boot_B, varargin)
% COMPARESTATEARRAYPOSTERIORS  Quantitative comparison of two saSPT D posteriors.
%
%   out = compareStateArrayPosteriors(boot_A, boot_B)
%   out = compareStateArrayPosteriors(boot_A, boot_B, ...
%             'SlowCutoff', 0.05, 'FastCutoff', 0.5, 'AlphaCI', 0.05)
%
%   Inputs
%     boot_A, boot_B : output structs from saSPT_bootstrap_CI.  Each must
%                      have fields .D, .post_marg_D_full, .post_marg_D_boot.
%                      The D grids must be identical.
%
%   Options
%     'SlowCutoff'  0.05 µm²/s  upper bound for the "slow/bound" fraction
%     'FastCutoff'  0.5  µm²/s  lower bound for the "fast/free" fraction
%     'AlphaCI'     0.05        percentile CI level (95% CI by default)
%     'Seed'        []          rng seed for the random pairing of
%                               bootstrap iterations across conditions
%
%   The CIs are computed by pairing bootstrap iterations across the two
%   conditions (random unmatched pairing) and re-evaluating each metric
%   per pair, so they reflect the joint variability rather than treating
%   the two conditions as independent.  Differences (e.g. slowFrac_diff)
%   are computed pairwise inside the loop, which gives a tighter CI than
%   subtracting independent CIs.
%
%   Output struct
%     D                       shared D grid (column)
%     post_A, post_B          full-pool point posterior for each condition
%     W1_log10D               Wasserstein-1 on log10(D), in decades
%     JS                      Jensen-Shannon divergence (nats)
%     JS_metric               sqrt(JS), the JS metric (a true distance)
%     slowFrac_A, slowFrac_B  occupancy below SlowCutoff
%     slowFrac_diff           slowFrac_A - slowFrac_B (paired)
%     fastFrac_A, fastFrac_B  occupancy above FastCutoff
%     fastFrac_diff           fastFrac_A - fastFrac_B (paired)
%     peakD_A, peakD_B        argmax D for each condition
%     peakShift_log10         log10(peakD_B / peakD_A) (paired)
%     opts                    the input options used
%
%   Each metric field is itself a struct with .point (full-pool point
%   estimate) and .CI = [lo hi] (percentile CI from the bootstrap pairing).
%
%   Recommendation: report W1 (single-number summary) plus the
%   bound/free-fraction differences (biological interpretation), all with
%   their CIs.

    p = inputParser;
    addParameter(p, 'SlowCutoff', 0.05, @(x) isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p, 'FastCutoff', 0.5,  @(x) isnumeric(x)&&isscalar(x)&&x>0);
    addParameter(p, 'AlphaCI',    0.05, @(x) isnumeric(x)&&isscalar(x)&&x>0&&x<1);
    addParameter(p, 'Seed',       [],   @(x) isempty(x)||(isnumeric(x)&&isscalar(x)));
    parse(p, varargin{:});
    o = p.Results;

    %% Validate inputs
    required = {'D','post_marg_D_full','post_marg_D_boot'};
    for i = 1:numel(required)
        f = required{i};
        if ~isfield(boot_A, f) || ~isfield(boot_B, f)
            error('compareStateArrayPosteriors: both inputs need field ''%s''.', f);
        end
    end
    if numel(boot_A.D) ~= numel(boot_B.D) || ...
       any(abs(boot_A.D(:) - boot_B.D(:)) > 1e-12 * max(1, abs(boot_A.D(:))))
        error('compareStateArrayPosteriors: D grids do not match between the two bootstraps.');
    end

    D    = boot_A.D(:);
    logD = log10(D);
    pA_full = normalize_(boot_A.post_marg_D_full);
    pB_full = normalize_(boot_B.post_marg_D_full);

    %% Pair bootstrap iterations for joint CI
    BA = size(boot_A.post_marg_D_boot, 2);
    BB = size(boot_B.post_marg_D_boot, 2);
    B  = min(BA, BB);
    if B < 10
        warning('compareStateArrayPosteriors:fewSamples', ...
            'Only %d paired bootstrap samples; CIs will be wide.', B);
    end
    if ~isempty(o.Seed), rng(o.Seed); end
    permA = randperm(BA, B);
    permB = randperm(BB, B);

    W1  = nan(B,1);
    JS  = nan(B,1);
    sFA = nan(B,1); sFB = nan(B,1); sFD = nan(B,1);
    fFA = nan(B,1); fFB = nan(B,1); fFD = nan(B,1);
    pkA = nan(B,1); pkB = nan(B,1); pkS = nan(B,1);

    slowMask = D < o.SlowCutoff;
    fastMask = D > o.FastCutoff;

    for b = 1:B
        pA = normalize_(boot_A.post_marg_D_boot(:, permA(b)));
        pB = normalize_(boot_B.post_marg_D_boot(:, permB(b)));
        if any(isnan(pA)) || any(isnan(pB)) || sum(pA) <= 0 || sum(pB) <= 0
            continue;
        end
        W1(b)  = wasserstein1(pA, pB, logD);
        JS(b)  = jensen_shannon(pA, pB);
        sFA(b) = sum(pA(slowMask));  sFB(b) = sum(pB(slowMask));  sFD(b) = sFA(b) - sFB(b);
        fFA(b) = sum(pA(fastMask));  fFB(b) = sum(pB(fastMask));  fFD(b) = fFA(b) - fFB(b);
        [~, iA] = max(pA); pkA(b) = D(iA);
        [~, iB] = max(pB); pkB(b) = D(iB);
        pkS(b) = log10(pkB(b) / pkA(b));
    end

    %% Point estimates from the full-pool posterior
    W1_pt  = wasserstein1(pA_full, pB_full, logD);
    JS_pt  = jensen_shannon(pA_full, pB_full);
    sFA_pt = sum(pA_full(slowMask));  sFB_pt = sum(pB_full(slowMask));
    fFA_pt = sum(pA_full(fastMask));  fFB_pt = sum(pB_full(fastMask));
    [~, iA] = max(pA_full); pkA_pt = D(iA);
    [~, iB] = max(pB_full); pkB_pt = D(iB);

    qLo = 100 * o.AlphaCI/2;
    qHi = 100 * (1 - o.AlphaCI/2);
    ciFn = @(v) prctile(v(~isnan(v)), [qLo qHi]);

    JS_clip = max(JS, 0);   % numerical safety for sqrt

    out = struct( ...
        'D',               D, ...
        'post_A',          pA_full, ...
        'post_B',          pB_full, ...
        'W1_log10D',       metric_(W1_pt,                    ciFn(W1)), ...
        'JS',              metric_(JS_pt,                    ciFn(JS)), ...
        'JS_metric',       metric_(sqrt(max(JS_pt,0)),       sqrt(ciFn(JS_clip))), ...
        'slowFrac_A',      metric_(sFA_pt,                   ciFn(sFA)), ...
        'slowFrac_B',      metric_(sFB_pt,                   ciFn(sFB)), ...
        'slowFrac_diff',   metric_(sFA_pt - sFB_pt,          ciFn(sFD)), ...
        'fastFrac_A',      metric_(fFA_pt,                   ciFn(fFA)), ...
        'fastFrac_B',      metric_(fFB_pt,                   ciFn(fFB)), ...
        'fastFrac_diff',   metric_(fFA_pt - fFB_pt,          ciFn(fFD)), ...
        'peakD_A',         metric_(pkA_pt,                   ciFn(pkA)), ...
        'peakD_B',         metric_(pkB_pt,                   ciFn(pkB)), ...
        'peakShift_log10', metric_(log10(pkB_pt/pkA_pt),     ciFn(pkS)), ...
        'opts',            o);
end

% -------------------------------------------------------------------------

function s = metric_(point, ci)
    s = struct('point', point, 'CI', ci);
end

function p = normalize_(p)
    p = p(:);
    s = sum(p);
    if s > 0, p = p / s; end
end

function w = wasserstein1(p, q, x)
    % 1D Wasserstein-1 distance between two discrete densities sharing
    % the support x: integral of |F_p - F_q| over x.
    cdfP = cumsum(p);
    cdfQ = cumsum(q);
    w = trapz(x, abs(cdfP - cdfQ));
end

function js = jensen_shannon(p, q)
    m  = 0.5 * (p + q);
    js = 0.5 * kl_safe(p, m) + 0.5 * kl_safe(q, m);
end

function k = kl_safe(p, q)
    % KL(p||q), with 0*log(0) = 0 and a tiny floor on q to avoid -inf
    % when bins of q are zero where p is also zero (no contribution).
    mask = p > 0;
    if ~any(mask)
        k = 0; return;
    end
    pp = p(mask);
    qq = q(mask) + 1e-300;
    k  = sum(pp .* (log(pp) - log(qq)));
end

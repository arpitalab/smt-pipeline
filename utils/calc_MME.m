function results = calc_MME(tracks, varargin)
%CALC_MME  Mean Maximal Excursion (MME) analysis for 2-D single-particle trajectories.
%
%   Implements Tejedor et al., Biophys. J. 98:1364-1372 (2010).
%
%   SYNTAX
%       results = calc_MME(tracks)
%       results = calc_MME(tracks, 'Param', value, ...)
%
%   INPUT
%       tracks    Cell array. Each cell is an Nx3 matrix [x, y, frame].
%                 x and y are positions (any consistent unit, e.g. um).
%                 frame is the integer frame number; gaps are allowed.
%
%   PARAMETERS
%       'MaxLag'            Maximum lag in frames. Default: min(floor(max_span/2), 100),
%                           at least 10.
%       'MinN'              Min tracks / TA windows per lag to include in fit.
%                           Default: 5.
%       'FitFrac'           [lo hi] fraction of the valid-lag range to use for
%                           power-law fitting. Default: [0.1 0.9].
%       'r0Scale'           Growing-sphere scale: r0 = r0Scale*mean_r(lag=1).
%                           Default: 2.
%       'Mode'              'auto' | 'ensemble' | 'timeaveraged' | 'both'
%                           'auto': uses 'timeaveraged' if mean track span >= 100
%                           frames, 'both' otherwise.
%                           Ensemble moments are ALWAYS computed (needed for
%                           moment ratios and classification). Mode controls
%                           whether TA moments are additionally computed.
%                           Default: 'auto'.
%       'Plot'              true | false. Default: true.
%
%   NOISE CORRECTION PARAMETERS
%       Real SMT data contains two sources of additive bias in MSD/MME2:
%         (1) Localization error (sigma_loc): adds +4*sigma^2 in 2D
%         (2) Motion blur (exposure time tau): subtracts a correction term
%       Both appear as an additive offset B in: MSD_meas(lag) = K*lag^alpha + B
%       where B>0 if localization dominates, B<0 if motion blur dominates.
%
%       'LocalizationError'  Known localization precision sigma (same units as
%                            x,y positions). Subtracts 4*sigma^2 from MSD/MME2
%                            before fitting alpha. Default: 0 (no correction).
%       'FitOffset'          true | false. Estimates B as a free parameter by
%                            scanning B values and minimizing log-log residuals.
%                            Ignored if LocalizationError > 0. Default: false.
%
%   OUTPUT  (struct)
%       .mode, .N_tracks, .mean_length, .lags
%       .ens              Ensemble results (always present):
%           .MSD, .MME2           Raw (uncorrected) moments
%           .MSD_fit, .MME2_fit   Noise-corrected versions used for fitting
%           .reg_ratio, .mme_ratio, .Neff
%           .alpha_MSD, .alpha_MME
%       .ta               Time-averaged results (when mode includes TA):
%           .MSD, .MME2, .MSD_fit, .MME2_fit, .Neff
%           .alpha_MSD, .alpha_MME
%       .alpha_MSD        Primary alpha (TA if available, else ensemble)
%       .alpha_MME        Primary alpha' from MME2
%       .alpha_corrected  FBM-corrected: (alpha_MME - 0.156)/0.849  [Eq.10]
%       .noise_correction
%           .mode         'none' | 'sigma_provided' | 'offset_fitted'
%           .B_MSD_ens    Offset applied to ensemble MSD
%           .B_MME2_ens   Offset applied to ensemble MME2
%           .B_MSD_ta     Offset applied to TA MSD (NaN if no TA)
%           .B_MME2_ta    Offset applied to TA MME2 (NaN if no TA)
%           .sigma_loc    Effective sigma_loc = sqrt(max(B_MSD,0)/4) in 2D
%       .growing_sphere   .prob, .lags, .slope, .df, .r0
%       .moment_ratios    .obs_reg, .obs_mme, .bm_reg, .bm_mme,
%                         .ctrw_reg, .ctrw_mme, .fbm_mme, .alpha_used
%       .classification   .process, .confidence, .notes
%
%   REFERENCE
%       Tejedor V et al. (2010) Biophys. J. 98:1364-1372.

%% Parse inputs
p = inputParser;
addRequired(p,  'tracks');
addParameter(p, 'MaxLag',            [],    @(x) isscalar(x) && x > 0);
addParameter(p, 'MinN',              5,     @(x) isscalar(x) && x > 0);
addParameter(p, 'FitFrac',           [0.1 0.9], @(x) numel(x)==2 && all(x>0) && x(1)<x(2));
addParameter(p, 'r0Scale',           2,     @isscalar);
addParameter(p, 'Mode',              'auto',@(x) ismember(lower(x), ...
    {'auto','ensemble','timeaveraged','both'}));
addParameter(p, 'Plot',              true,  @(x) islogical(x) || ismember(x,[0 1]));
addParameter(p, 'LocalizationError', 0,     @(x) isscalar(x) && x >= 0);
addParameter(p, 'FitOffset',         false, @(x) islogical(x) || ismember(x,[0 1]));
addParameter(p, 'Bootstrap',         0,     @(x) isscalar(x) && x >= 0 && x == round(x));
parse(p, tracks, varargin{:});
opt = p.Results;

%% Validate tracks
nPts = cellfun('size', tracks, 1);
bad  = nPts < 2;
if any(bad)
    warning('calc_MME: %d track(s) with fewer than 2 points removed.', sum(bad));
    tracks = tracks(~bad);
end
nTracks = numel(tracks);
if nTracks == 0
    error('calc_MME: No valid tracks (need >= 2 data points per track).');
end

%% Basic statistics and mode
spans    = cellfun(@(T) T(end,3) - T(1,3), tracks);
meanSpan = mean(spans);

switch lower(opt.Mode)
    case 'auto'
        if meanSpan >= 100, mode = 'timeaveraged'; else, mode = 'both'; end
    otherwise
        mode = lower(opt.Mode);
end

if isempty(opt.MaxLag)
    maxLag = max(10, min(floor(max(spans)/2), 100));
else
    maxLag = round(opt.MaxLag);
end
lags = (1:maxLag)';

doTA = ismember(mode, {'timeaveraged','both'});

%% Pre-process
proc = preprocess(tracks, nTracks);

%% Compute raw moments (always ensemble; TA only when requested)
[MSD_e, MME2_e, r4_e, r4max_e, Neff_e] = ensembleMoments(proc, lags);

if doTA
    [MSD_t, MME2_t, Neff_t] = taMoments(proc, lags);
end

%% Noise correction
% Determine offset B for each curve. B is subtracted before fitting alpha.
% Raw moment arrays are preserved unchanged in results for plotting.
if opt.LocalizationError > 0
    B = 4 * opt.LocalizationError^2;    % 2D: 2 axes x 2*sigma^2
    B_MSD_e  = B;   B_MME2_e = B;
    B_MSD_t  = B;   B_MME2_t = B;
    nc_mode  = 'sigma_provided';
    sigma_nc = opt.LocalizationError;

elseif opt.FitOffset
    B_MSD_e  = estimateOffset(lags, MSD_e,  Neff_e, opt.MinN, opt.FitFrac);
    B_MME2_e = estimateOffset(lags, MME2_e, Neff_e, opt.MinN, opt.FitFrac);
    if doTA
        B_MSD_t  = estimateOffset(lags, MSD_t,  Neff_t, opt.MinN, opt.FitFrac);
        B_MME2_t = estimateOffset(lags, MME2_t, Neff_t, opt.MinN, opt.FitFrac);
    else
        B_MSD_t = NaN; B_MME2_t = NaN;
    end
    nc_mode  = 'offset_fitted';
    sigma_nc = sqrt(max(B_MSD_e, 0) / 4);   % B = 4*sigma^2 in 2D

else
    B_MSD_e = 0;  B_MME2_e = 0;
    B_MSD_t = 0;  B_MME2_t = 0;
    nc_mode  = 'none';
    sigma_nc = NaN;
end

% Apply offsets: subtract B, replace non-positive results with NaN
MSD_e_fit  = floorNaN(MSD_e,  B_MSD_e);
MME2_e_fit = floorNaN(MME2_e, B_MME2_e);
if doTA
    MSD_t_fit  = floorNaN(MSD_t,  B_MSD_t);
    MME2_t_fit = floorNaN(MME2_t, B_MME2_t);
end

%% Fit power laws using corrected data
rr_e = r4_e    ./ (MSD_e.^2);
mr_e = r4max_e ./ (MME2_e.^2);

[aM_e, ~] = fitPL(lags, MSD_e_fit,  Neff_e, opt.MinN, opt.FitFrac);
[aP_e, ~] = fitPL(lags, MME2_e_fit, Neff_e, opt.MinN, opt.FitFrac);

results.ens = struct('MSD', MSD_e, 'MME2', MME2_e, ...
                     'MSD_fit', MSD_e_fit, 'MME2_fit', MME2_e_fit, ...
                     'reg_ratio', rr_e, 'mme_ratio', mr_e, ...
                     'Neff', Neff_e, 'alpha_MSD', aM_e, 'alpha_MME', aP_e);

if doTA
    [aM_t, ~] = fitPL(lags, MSD_t_fit,  Neff_t, opt.MinN, opt.FitFrac);
    [aP_t, ~] = fitPL(lags, MME2_t_fit, Neff_t, opt.MinN, opt.FitFrac);
    results.ta = struct('MSD', MSD_t, 'MME2', MME2_t, ...
                        'MSD_fit', MSD_t_fit, 'MME2_fit', MME2_t_fit, ...
                        'Neff', Neff_t, 'alpha_MSD', aM_t, 'alpha_MME', aP_t);
end

%% Fill results header
results.mode        = mode;
results.N_tracks    = nTracks;
results.mean_length = meanSpan;
results.lags        = lags;

results.noise_correction = struct( ...
    'mode',       nc_mode, ...
    'B_MSD_ens',  B_MSD_e, ...
    'B_MME2_ens', B_MME2_e, ...
    'B_MSD_ta',   B_MSD_t, ...
    'B_MME2_ta',  B_MME2_t, ...
    'sigma_loc',  sigma_nc);

%% Primary alpha (TA preferred for long tracks)
if doTA
    alpha_MSD = results.ta.alpha_MSD;
    alpha_MME = results.ta.alpha_MME;
else
    alpha_MSD = results.ens.alpha_MSD;
    alpha_MME = results.ens.alpha_MME;
end
results.alpha_MSD       = alpha_MSD;
results.alpha_MME       = alpha_MME;
results.alpha_corrected = (alpha_MME - 0.156) / 0.849;

%% Growing sphere (uses ensemble at exact lags; no noise correction needed
%  because r_max is taken from raw positions, not displacements)
r0 = opt.r0Scale * meanR1(proc, lags);
results.growing_sphere = growingSphere(proc, lags, alpha_MSD, r0, opt.MinN);

%% Classification (always uses ensemble moment ratios)
[results.moment_ratios, results.classification] = classifyProcess( ...
    results.ens.reg_ratio, results.ens.mme_ratio, results.ens.Neff, ...
    alpha_MSD, alpha_MME, results.growing_sphere, opt.MinN);

%% Bootstrap CI (optional)
% Resamples tracks with replacement and recomputes all fitted quantities
% to obtain 95% CI (2.5–97.5 percentile) on alpha_MSD, alpha_MME,
% moment ratios, noise offsets, growing-sphere slope, and pointwise
% CI bands on MSD/MME2 curves.
if opt.Bootstrap > 0
    fprintf('calc_MME: bootstrapping (%d iterations) ', opt.Bootstrap);
    results.bootstrap = bootstrapCI(proc, lags, doTA, opt, round(opt.Bootstrap));
    fprintf(' done.\n');
else
    results.bootstrap = struct('nBoot', 0);
end

%% Plots
if opt.Plot
    makePlots(results, mode);
end

end  % main


%% =========================================================================
function proc = preprocess(tracks, nTracks)
proc(nTracks) = struct('relFrames',[],'dists',[],'x',[],'y',[],'frames',[]);
for i = 1:nTracks
    T = tracks{i};
    proc(i).frames    = T(:,3);
    proc(i).x         = T(:,1);
    proc(i).y         = T(:,2);
    proc(i).relFrames = T(:,3) - T(1,3);
    proc(i).dists     = sqrt((T(:,1)-T(1,1)).^2 + (T(:,2)-T(1,2)).^2);
end
end


%% =========================================================================
%%  ENSEMBLE MOMENTS  (Eqs. 17-18)
%% =========================================================================
function [MSD, MME2, r4, r4max, Neff] = ensembleMoments(proc, lags)
nT   = numel(proc);
nL   = numel(lags);
MSD  = nan(nL,1);  MME2  = nan(nL,1);
r4   = nan(nL,1);  r4max = nan(nL,1);
Neff = zeros(nL,1);

for li = 1:nL
    t  = lags(li);
    rv = nan(nT,1);
    rm = nan(nT,1);
    for i = 1:nT
        rf = proc(i).relFrames;
        d  = proc(i).dists;
        k  = find(rf == t, 1);
        if isempty(k), continue; end
        rv(i) = d(k);
        rm(i) = max(d(rf <= t));
    end
    ok       = ~isnan(rv);
    Neff(li) = sum(ok);
    if Neff(li) < 2, continue; end
    r  = rv(ok);  rmv = rm(ok);
    MSD(li)   = mean(r.^2);
    MME2(li)  = mean(rmv.^2);
    r4(li)    = mean(r.^4);
    r4max(li) = mean(rmv.^4);
end
end


%% =========================================================================
%%  TIME-AVERAGED MOMENTS  (Eqs. 2, 20)
%% =========================================================================
function [taMSD, taMME2, Neff] = taMoments(proc, lags)
nL      = numel(lags);
sumMSD  = zeros(nL,1);
sumMME2 = zeros(nL,1);
nWin    = zeros(nL,1);

for i = 1:numel(proc)
    fr = proc(i).frames;  x = proc(i).x;  y = proc(i).y;
    nP = numel(fr);
    if nP < 2, continue; end
    f0 = fr(1);  fN = fr(end);
    if fN - f0 < lags(1), continue; end

    lut = zeros(fN - f0 + 1, 1, 'int32');
    for k = 1:nP
        lut(fr(k) - f0 + 1) = k;
    end

    for ji = 1:nP
        fji   = fr(ji);
        maxSq = 0;
        prevE = ji;

        for li = 1:nL
            fEnd = fji + lags(li);
            if fEnd > fN, break; end
            ei = lut(fEnd - f0 + 1);
            if ei == 0, continue; end

            for k = (prevE+1):ei
                d2 = (x(k)-x(ji))^2 + (y(k)-y(ji))^2;
                if d2 > maxSq, maxSq = d2; end
            end
            prevE = ei;

            sumMSD(li)  = sumMSD(li)  + (x(ei)-x(ji))^2 + (y(ei)-y(ji))^2;
            sumMME2(li) = sumMME2(li) + maxSq;
            nWin(li)    = nWin(li) + 1;
        end
    end
end

Neff   = nWin;
taMSD  = nan(nL,1);  taMME2 = nan(nL,1);
ok = nWin > 0;
taMSD(ok)  = sumMSD(ok)  ./ nWin(ok);
taMME2(ok) = sumMME2(ok) ./ nWin(ok);
end


%% =========================================================================
%%  NOISE CORRECTION HELPERS
%% =========================================================================

% Estimate additive offset B by scanning over B values and minimizing
% log-log residuals of (vals - B) vs lags.  Uses only the fit-range lags.
function B = estimateOffset(lags, vals, Neff, minN, fitFrac)
valid = Neff >= minN & isfinite(vals) & vals > 0;
idx   = find(valid);
if numel(idx) < 4, B = 0; return; end

nv  = numel(idx);
lo  = idx(max(1,   round(fitFrac(1)*nv)));
hi  = idx(min(nv,  round(fitFrac(2)*nv)));
mask = false(numel(lags), 1);
mask(lo:hi) = true;
mask = mask & valid;

lv = lags(mask);
vv = vals(mask);
if numel(lv) < 3, B = 0; return; end

% Search range: B in (-(max-min), 0.9*min)
% Negative B: motion blur dominated (MSD at lag=1 lower than expected)
% Positive B: localization dominated (MSD floor at short lags)
B_min = -(max(vv) - min(vv));
B_max = 0.90 * min(vv);

obj = @(b) offsetObjective(lv, vv, b);
B   = fminbnd(obj, B_min, B_max);

% If B is negligible relative to signal, snap to zero
if abs(B) < 0.01 * mean(vv), B = 0; end
end

function res = offsetObjective(lv, vv, B)
pos = vv - B;
if any(pos <= 0), res = 1e10; return; end
p    = polyfit(log(lv), log(pos), 1);
res  = mean((log(pos) - polyval(p, log(lv))).^2);
end

% Subtract offset and replace non-positive entries with NaN
function out = floorNaN(vals, B)
if B == 0
    out = vals;
    return;
end
out = vals - B;
out(out <= 0) = NaN;
end


%% =========================================================================
%%  POWER-LAW FIT  (log-log linear regression)
%% =========================================================================
function [alpha, logK] = fitPL(lags, vals, Neff, minN, fitFrac)
valid = Neff >= minN & isfinite(vals) & vals > 0;
idx   = find(valid);
if numel(idx) < 3, alpha = NaN; logK = NaN; return; end

nv   = numel(idx);
lo   = idx(max(1,   round(fitFrac(1)*nv)));
hi   = idx(min(nv,  round(fitFrac(2)*nv)));
mask = false(numel(lags), 1);
mask(lo:hi) = true;
mask = mask & valid;
if sum(mask) < 2, alpha = NaN; logK = NaN; return; end

p     = polyfit(log(lags(mask)), log(vals(mask)), 1);
alpha = p(1);
logK  = p(2);
end


%% =========================================================================
%%  MEAN DISTANCE AT FIRST LAG
%% =========================================================================
function r1 = meanR1(proc, lags)
t  = lags(1);
rv = [];
for i = 1:numel(proc)
    k = find(proc(i).relFrames == t, 1);
    if ~isempty(k), rv(end+1) = proc(i).dists(k); end %#ok<AGROW>
end
r1 = max(mean(rv(isfinite(rv))), 1e-10);
end


%% =========================================================================
%%  GROWING SPHERE  (Eq. 19)
%% =========================================================================
function gs = growingSphere(proc, lags, alpha, r0, minN)
nL   = numel(lags);
prob = nan(nL,1);
Neff = zeros(nL,1);

if isnan(alpha), alpha_gs = 1; else, alpha_gs = alpha; end

for li = 1:nL
    t    = lags(li);
    rsph = r0 * t^(alpha_gs/2);
    nIn  = 0;  nTot = 0;
    for i = 1:numel(proc)
        k = find(proc(i).relFrames == t, 1);
        if isempty(k), continue; end
        nTot = nTot + 1;
        if proc(i).dists(k) <= rsph, nIn = nIn + 1; end
    end
    Neff(li) = nTot;
    if nTot >= minN, prob(li) = nIn / nTot; end
end

valid = isfinite(prob) & prob > 0 & Neff >= minN;
if sum(valid) >= 3
    pp    = polyfit(log(lags(valid)), log(prob(valid)), 1);
    slope = pp(1);
    df    = 2 - 2*slope/alpha_gs;
else
    slope = NaN;  df = NaN;
end

gs = struct('prob',prob,'lags',lags,'slope',slope,'df',df,'r0',r0);
end


%% =========================================================================
%%  PROCESS CLASSIFICATION  (Table 2)
%% =========================================================================
function [mr, cls] = classifyProcess(reg_ratio, mme_ratio, Neff, ...
                                      alpha_MSD, alpha_MME, gs, minN)
bm_reg = 2.00;
bm_mme = 1.49;

if ~isnan(alpha_MSD)
    a     = max(0.05, alpha_MSD);
    scale = 2 * gamma(a+1)^2 / gamma(2*a+1);
    ctrw_reg = scale * bm_reg;
    ctrw_mme = scale * bm_mme;
    fbm_mme  = 1.055 * (a/2)^1.42 + 1.10;
else
    [ctrw_reg, ctrw_mme, fbm_mme] = deal(NaN);
end

[obs_reg, obs_mme] = deal(NaN);
if ~isempty(reg_ratio) && ~isempty(Neff)
    nL   = numel(reg_ratio);
    tail = max(1, round(0.70*nL)):nL;
    ok   = tail(Neff(tail) >= minN);
    if numel(ok) >= 3
        obs_reg = nanmedian(reg_ratio(ok));
        obs_mme = nanmedian(mme_ratio(ok));
    end
end

mr = struct('obs_reg',obs_reg,'obs_mme',obs_mme,'bm_reg',bm_reg,'bm_mme',bm_mme,...
            'ctrw_reg',ctrw_reg,'ctrw_mme',ctrw_mme,'fbm_mme',fbm_mme,...
            'alpha_used',alpha_MSD);

tol   = 0.25;   % tolerance for "consistent with BM reference value"
notes = {};
if isnan(alpha_MSD)
    cls = mkCls('undetermined','low',{'Could not fit alpha from MSD.'});
    return
end
if alpha_MSD > 0.85
    % Mostly normal diffusion, but still check moment ratios for consistency.
    % Near alpha=1 the CTRW and FBM ratios converge to BM values, so genuine
    % distinguishability is limited — we flag discrepancies rather than
    % reclassifying.
    bm_notes = {sprintf('alpha_MSD = %.2f, consistent with normal diffusion', alpha_MSD)};
    if ~isnan(obs_reg) && ~isnan(obs_mme)
        if obs_reg > bm_reg + tol && obs_mme > bm_mme + tol
            bm_notes{end+1} = sprintf(['WARNING: moment ratios elevated ' ...
                '(reg=%.2f, MME=%.2f) — may indicate weakly subdiffusive CTRW'], ...
                obs_reg, obs_mme);
            bm_notes{end+1} = sprintf(['Expected CTRW at this alpha: ' ...
                'reg=%.2f, MME=%.2f'], ctrw_reg, ctrw_mme);
        elseif abs(obs_reg - bm_reg) <= tol && obs_mme < bm_mme - tol
            bm_notes{end+1} = sprintf(['WARNING: MME ratio depressed ' ...
                '(MME=%.2f < 1.49) — may indicate weakly subdiffusive FBM'], ...
                obs_mme);
        else
            bm_notes{end+1} = sprintf('Moment ratios consistent with BM (reg=%.2f, MME=%.2f)', ...
                obs_reg, obs_mme);
        end
    end
    cls = mkCls('Brownian motion (BM)','high', bm_notes);
    return
end

is_fractal_gs = ~isnan(gs.slope) && gs.slope > 0.05;
if is_fractal_gs
    notes{end+1} = sprintf('Growing-sphere probability increases: slope = %.3f', gs.slope);
    if ~isnan(gs.df)
        notes{end+1} = sprintf('Estimated fractal dimension d_f = %.2f (d=2)', gs.df);
    end
    cls = mkCls('diffusion on fractal','medium',notes);
    return
end

if ~isnan(obs_reg)
    reg_above = obs_reg > bm_reg + tol;
    mme_above = obs_mme > bm_mme + tol;
    reg_near  = abs(obs_reg - bm_reg) <= tol;
    mme_below = obs_mme < bm_mme - tol;

    if reg_above && mme_above
        notes{end+1} = sprintf('reg = %.2f > 2.00, MME = %.2f > 1.49', obs_reg, obs_mme);
        notes{end+1} = sprintf('Expected CTRW: reg = %.2f, MME = %.2f', ctrw_reg, ctrw_mme);
        notes{end+1} = 'Growing-sphere is constant, excluding fractal diffusion.';
        cls = mkCls('continuous time random walk (CTRW)','high',notes);
    elseif reg_near && mme_below
        notes{end+1} = sprintf('reg = %.2f ≈ 2.00 (BM)', obs_reg);
        notes{end+1} = sprintf('MME = %.2f < 1.49 (BM), expected FBM = %.2f', obs_mme, fbm_mme);
        if ~isnan(alpha_MME) && alpha_MME > alpha_MSD + 0.03
            notes{end+1} = sprintf('alpha_MME = %.2f > alpha_MSD = %.2f (expected for FBM)', ...
                                    alpha_MME, alpha_MSD);
            notes{end+1} = sprintf('FBM-corrected alpha = %.2f', (alpha_MME-0.156)/0.849);
        end
        cls = mkCls('fractional Brownian motion (FBM)','high',notes);
    elseif obs_reg < bm_reg - tol
        notes{end+1} = sprintf('reg = %.2f < 2.00, consistent with fractal (possible confinement)', obs_reg);
        cls = mkCls('diffusion on fractal','medium',notes);
    else
        notes{end+1} = sprintf('reg = %.2f, MME = %.2f — ratios ambiguous', obs_reg, obs_mme);
        notes{end+1} = sprintf('BM: (%.2f,%.2f)  CTRW: (%.2f,%.2f)  FBM mme: %.2f', ...
                                bm_reg,bm_mme,ctrw_reg,ctrw_mme,fbm_mme);
        cls = mkCls('subdiffusion (ambiguous moment ratios)','low',notes);
    end
else
    if ~isnan(alpha_MME) && ~isnan(alpha_MSD)
        fbm_prime = 0.156 + 0.849 * alpha_MSD;
        if alpha_MME > alpha_MSD + 0.03 && abs(alpha_MME - fbm_prime) < 0.10
            notes{end+1} = sprintf('alpha_MME = %.2f consistent with FBM alpha'' = %.2f', ...
                                    alpha_MME, fbm_prime);
            notes{end+1} = sprintf('FBM-corrected alpha = %.2f', (alpha_MME-0.156)/0.849);
            cls = mkCls('fractional Brownian motion (FBM) [indicative]','medium',notes);
        else
            notes{end+1} = 'Mechanism undetermined without moment ratios.';
            notes{end+1} = 'Re-run with Mode=''ensemble'' or ''both''.';
            cls = mkCls('subdiffusion (mechanism undetermined)','low',notes);
        end
    else
        cls = mkCls('undetermined','low',{'Could not fit alpha.'});
    end
end
end

function cls = mkCls(process, confidence, notes)
cls = struct('process',process,'confidence',confidence,'notes',{notes});
end


%% =========================================================================
%%  BOOTSTRAP CI
%% =========================================================================
%  Resamples nT tracks with replacement on each iteration.
%  This preserves within-track temporal correlations, which is important
%  for the time-averaged computation where windows from the same track are
%  not independent.
%  NOTE: for TA mode with long tracks, each iteration calls taMoments
%  (the slow nested loop).  Keep nBoot <= 500 for reasonable runtimes,
%  or use Mode='ensemble' for faster bootstrapping.
function boot = bootstrapCI(proc, lags, doTA, opt, nBoot)
nT   = numel(proc);
nL   = numel(lags);
tail = max(1, round(0.70*nL));   % long-time plateau start index

% Pre-allocate scalar distributions
alpha_MSD_b  = nan(nBoot,1);
alpha_MME_b  = nan(nBoot,1);
obs_reg_b    = nan(nBoot,1);
obs_mme_b    = nan(nBoot,1);
B_MSD_b      = nan(nBoot,1);
B_MME2_b     = nan(nBoot,1);
gs_slope_b   = nan(nBoot,1);

% Pre-allocate per-lag distributions for pointwise CI bands
MSD_ens_b    = nan(nBoot,nL);
MME2_ens_b   = nan(nBoot,nL);
MSD_ta_b     = nan(nBoot,nL);
MME2_ta_b    = nan(nBoot,nL);

printEvery = max(1, round(nBoot/10));

for b = 1:nBoot
    % ── Resample tracks ──────────────────────────────────────────────────
    idx    = randi(nT, nT, 1);
    proc_b = proc(idx);

    % ── Ensemble moments ──────────────────────────────────────────────────
    [MSD_e, MME2_e, r4_e, r4max_e, Neff_e] = ensembleMoments(proc_b, lags);
    MSD_ens_b(b,:)  = MSD_e';
    MME2_ens_b(b,:) = MME2_e';

    % ── TA moments ────────────────────────────────────────────────────────
    if doTA
        [MSD_t, MME2_t, Neff_t] = taMoments(proc_b, lags);
        MSD_ta_b(b,:)  = MSD_t';
        MME2_ta_b(b,:) = MME2_t';
    end

    % ── Noise correction ─────────────────────────────────────────────────
    if opt.LocalizationError > 0
        Bval         = 4 * opt.LocalizationError^2;
        MSD_e_fit    = floorNaN(MSD_e,  Bval);
        MME2_e_fit   = floorNaN(MME2_e, Bval);
        B_MSD_b(b)   = Bval;
        B_MME2_b(b)  = Bval;
        if doTA
            MSD_t_fit  = floorNaN(MSD_t,  Bval);
            MME2_t_fit = floorNaN(MME2_t, Bval);
        end
    elseif opt.FitOffset
        Bm           = estimateOffset(lags, MSD_e,  Neff_e, opt.MinN, opt.FitFrac);
        Bm2          = estimateOffset(lags, MME2_e, Neff_e, opt.MinN, opt.FitFrac);
        MSD_e_fit    = floorNaN(MSD_e,  Bm);
        MME2_e_fit   = floorNaN(MME2_e, Bm2);
        B_MSD_b(b)   = Bm;
        B_MME2_b(b)  = Bm2;
        if doTA
            Bmt         = estimateOffset(lags, MSD_t,  Neff_t, opt.MinN, opt.FitFrac);
            Bmt2        = estimateOffset(lags, MME2_t, Neff_t, opt.MinN, opt.FitFrac);
            MSD_t_fit   = floorNaN(MSD_t,  Bmt);
            MME2_t_fit  = floorNaN(MME2_t, Bmt2);
        end
    else
        MSD_e_fit   = MSD_e;
        MME2_e_fit  = MME2_e;
        if doTA
            MSD_t_fit  = MSD_t;
            MME2_t_fit = MME2_t;
        end
    end

    % ── Fit alpha ─────────────────────────────────────────────────────────
    if doTA
        [aM, ~] = fitPL(lags, MSD_t_fit,  Neff_t, opt.MinN, opt.FitFrac);
        [aP, ~] = fitPL(lags, MME2_t_fit, Neff_t, opt.MinN, opt.FitFrac);
    else
        [aM, ~] = fitPL(lags, MSD_e_fit,  Neff_e, opt.MinN, opt.FitFrac);
        [aP, ~] = fitPL(lags, MME2_e_fit, Neff_e, opt.MinN, opt.FitFrac);
    end
    alpha_MSD_b(b) = aM;
    alpha_MME_b(b) = aP;

    % ── Moment ratios (ensemble, long-time plateau) ───────────────────────
    rr_b = r4_e ./ (MSD_e.^2);
    mr_b = r4max_e ./ (MME2_e.^2);
    ok   = (tail:nL);
    ok   = ok(Neff_e(ok) >= opt.MinN);
    if numel(ok) >= 3
        obs_reg_b(b) = nanmedian(rr_b(ok));
        obs_mme_b(b) = nanmedian(mr_b(ok));
    end

    % ── Growing sphere slope ──────────────────────────────────────────────
    r0_b         = opt.r0Scale * meanR1(proc_b, lags);
    gs_b         = growingSphere(proc_b, lags, aM, r0_b, opt.MinN);
    gs_slope_b(b) = gs_b.slope;

    if mod(b, printEvery) == 0, fprintf('.'); end
end

% ── Assemble output ───────────────────────────────────────────────────────
ci95 = @(v) prctile(v(isfinite(v)), [2.5 97.5]);

boot.nBoot = nBoot;

% Full distributions (useful for histograms / diagnostics)
boot.alpha_MSD_dist  = alpha_MSD_b;
boot.alpha_MME_dist  = alpha_MME_b;
boot.alpha_corr_dist = (alpha_MME_b - 0.156) / 0.849;
boot.obs_reg_dist    = obs_reg_b;
boot.obs_mme_dist    = obs_mme_b;
boot.B_MSD_dist      = B_MSD_b;
boot.B_MME2_dist     = B_MME2_b;
boot.gs_slope_dist   = gs_slope_b;

% 95% CIs (rows: [lo, hi])
boot.alpha_MSD_ci    = ci95(alpha_MSD_b);
boot.alpha_MME_ci    = ci95(alpha_MME_b);
boot.alpha_corr_ci   = ci95(boot.alpha_corr_dist);
boot.obs_reg_ci      = ci95(obs_reg_b);
boot.obs_mme_ci      = ci95(obs_mme_b);
boot.B_MSD_ci        = ci95(B_MSD_b);
boot.B_MME2_ci       = ci95(B_MME2_b);
boot.gs_slope_ci     = ci95(gs_slope_b);

% Standard errors
boot.alpha_MSD_se    = nanstd(alpha_MSD_b);
boot.alpha_MME_se    = nanstd(alpha_MME_b);
boot.obs_reg_se      = nanstd(obs_reg_b);
boot.obs_mme_se      = nanstd(obs_mme_b);

% Pointwise 95% CI bands on MSD/MME2 curves [nLags x 2]: col1=lo, col2=hi
boot.MSD_ens_ci  = [prctile(MSD_ens_b,   2.5, 1)', prctile(MSD_ens_b,  97.5, 1)'];
boot.MME2_ens_ci = [prctile(MME2_ens_b,  2.5, 1)', prctile(MME2_ens_b, 97.5, 1)'];
if doTA && any(isfinite(MSD_ta_b(:)))
    boot.MSD_ta_ci  = [prctile(MSD_ta_b,  2.5,1)', prctile(MSD_ta_b,  97.5,1)'];
    boot.MME2_ta_ci = [prctile(MME2_ta_b, 2.5,1)', prctile(MME2_ta_b, 97.5,1)'];
else
    boot.MSD_ta_ci  = nan(nL,2);
    boot.MME2_ta_ci = nan(nL,2);
end
end


%% =========================================================================
%%  PLOTS
%% =========================================================================
function makePlots(results, mode)
lags = results.lags;
mr   = results.moment_ratios;
gs   = results.growing_sphere;
cls  = results.classification;
nc   = results.noise_correction;

fig = figure('Name','MME Analysis','NumberTitle','off','Position',[60 60 1360 840]);

noiseStr = '';
if ~strcmp(nc.mode,'none')
    noiseStr = sprintf('  |  noise corr: %s (sigma~%.3g)', nc.mode, nc.sigma_loc);
end
sgtitle(sprintf('MME Analysis  |  N=%d  |  mean span=%.0f fr  |  mode: %s%s', ...
    results.N_tracks, results.mean_length, mode, noiseStr), ...
    'FontWeight','bold','FontSize',10);

%-- Panel 1: MSD and MME2 (log-log) ----------------------------------------
subplot(2,3,1); hold on; box on; grid on;

hasBoot = isfield(results,'bootstrap') && results.bootstrap.nBoot > 0;

% Bootstrap CI bands (drawn first so data sits on top)
if hasBoot
    boot = results.bootstrap;
    % Ensemble MSD CI band
    msd_lo = boot.MSD_ens_ci(:,1);
    msd_hi = boot.MSD_ens_ci(:,2);
    ok = isfinite(msd_lo) & isfinite(msd_hi) & msd_lo > 0 & msd_hi > 0;
    if any(ok)
        lv = lags(ok); ll = msd_lo(ok); lh = msd_hi(ok);
        fill([lv; flipud(lv)], [ll; flipud(lh)], [0 0.4 0.8], ...
             'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    end
    % Ensemble MME2 CI band
    mme_lo = boot.MME2_ens_ci(:,1);
    mme_hi = boot.MME2_ens_ci(:,2);
    ok = isfinite(mme_lo) & isfinite(mme_hi) & mme_lo > 0 & mme_hi > 0;
    if any(ok)
        lv = lags(ok); ll = mme_lo(ok); lh = mme_hi(ok);
        fill([lv; flipud(lv)], [ll; flipud(lh)], [0.8 0.2 0.2], ...
             'FaceAlpha',0.15,'EdgeColor','none','HandleVisibility','off');
    end
    % TA CI bands if present
    if isfield(boot,'MSD_ta_ci')
        msd_lo = boot.MSD_ta_ci(:,1); msd_hi = boot.MSD_ta_ci(:,2);
        ok = isfinite(msd_lo) & isfinite(msd_hi) & msd_lo > 0 & msd_hi > 0;
        if any(ok)
            lv = lags(ok);
            fill([lv; flipud(lv)], [msd_lo(ok); flipud(msd_hi(ok))], [0 0.4 0.8], ...
                 'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
        end
        mme_lo = boot.MME2_ta_ci(:,1); mme_hi = boot.MME2_ta_ci(:,2);
        ok = isfinite(mme_lo) & isfinite(mme_hi) & mme_lo > 0 & mme_hi > 0;
        if any(ok)
            lv = lags(ok);
            fill([lv; flipud(lv)], [mme_lo(ok); flipud(mme_hi(ok))], [0.8 0.2 0.2], ...
                 'FaceAlpha',0.10,'EdgeColor','none','HandleVisibility','off');
        end
    end
end

% Raw data
if isfield(results,'ens')
    loglog(lags, results.ens.MSD,  'bo','MarkerSize',4,'MarkerFaceColor','b', ...
           'DisplayName','MSD (ens, raw)');
    loglog(lags, results.ens.MME2, 'rs','MarkerSize',4,'MarkerFaceColor','r', ...
           'DisplayName','MME_2 (ens, raw)');
    % Fit lines use corrected data
    addFitLine(lags, results.ens.MSD_fit,  results.ens.alpha_MSD, ...
               nc.B_MSD_ens,  'b--', sprintf('\\alpha_{MSD}=%.2f (ens)',  results.ens.alpha_MSD));
    addFitLine(lags, results.ens.MME2_fit, results.ens.alpha_MME, ...
               nc.B_MME2_ens, 'r--', sprintf('\\alpha''_{MME}=%.2f (ens)', results.ens.alpha_MME));
end
if isfield(results,'ta')
    loglog(lags, results.ta.MSD,  'b^','MarkerSize',4,'DisplayName','MSD (TA, raw)');
    loglog(lags, results.ta.MME2, 'r^','MarkerSize',4,'DisplayName','MME_2 (TA, raw)');
    addFitLine(lags, results.ta.MSD_fit,  results.ta.alpha_MSD, ...
               nc.B_MSD_ta,  'b:', sprintf('\\alpha_{MSD}=%.2f (TA)',  results.ta.alpha_MSD));
    addFitLine(lags, results.ta.MME2_fit, results.ta.alpha_MME, ...
               nc.B_MME2_ta, 'r:', sprintf('\\alpha''_{MME}=%.2f (TA)', results.ta.alpha_MME));
end
% Annotate noise floor
if ~strcmp(nc.mode,'none') && nc.B_MSD_ens ~= 0
    yline(abs(nc.B_MSD_ens),'k:','LineWidth',1, ...
          'DisplayName',sprintf('|B|=%.3g (%s)', abs(nc.B_MSD_ens), ...
          ternary(nc.B_MSD_ens>0,'loc. floor','blur dip')));
end
xlabel('Lag (frames)'); ylabel('\langler^2\rangle');
title('MSD and 2^{nd} MME Moment');
legend('Location','northwest','FontSize',7);

%-- Panel 2: Regular moment ratio ------------------------------------------
subplot(2,3,2); hold on; box on; grid on;
if isfield(results,'ens')
    plot(lags, results.ens.reg_ratio, 'b.','MarkerSize',8);
end
yline(mr.bm_reg,'k-','LineWidth',1.5);
if ~isnan(mr.ctrw_reg), yline(mr.ctrw_reg,'m--','LineWidth',1.5); end
if ~isnan(mr.obs_reg),  yline(mr.obs_reg, 'r-', 'LineWidth',1.5); end
% Bootstrap CI for obs_reg: horizontal shaded band
if hasBoot && ~any(isnan(boot.obs_reg_ci))
    xl2 = xlim; if xl2(1)==0 && xl2(2)==1, xl2=[lags(1) lags(end)]; end
    fill([xl2(1) xl2(2) xl2(2) xl2(1)], ...
         [boot.obs_reg_ci(1) boot.obs_reg_ci(1) boot.obs_reg_ci(2) boot.obs_reg_ci(2)], ...
         [1 0 0],'FaceAlpha',0.12,'EdgeColor','none');
end
xl = xlim;
text(xl(2)*0.55, mr.bm_reg+0.12,'BM (2.00)','Color','k','FontSize',8);
if ~isnan(mr.ctrw_reg)
    text(xl(2)*0.55,mr.ctrw_reg+0.12,sprintf('CTRW (%.2f)',mr.ctrw_reg),'Color','m','FontSize',8);
end
if ~isnan(mr.obs_reg)
    text(xl(2)*0.55,mr.obs_reg+0.12,sprintf('Obs (%.2f)',mr.obs_reg),'Color','r','FontSize',8);
end
ylim([0 max(6, (mr.ctrw_reg+1)*~isnan(mr.ctrw_reg) + 6*isnan(mr.ctrw_reg))]);
xlabel('Lag (frames)'); ylabel('\langler^4\rangle/\langler^2\rangle^2');
title('Regular Moment Ratio');

%-- Panel 3: MME moment ratio -----------------------------------------------
subplot(2,3,3); hold on; box on; grid on;
if isfield(results,'ens')
    plot(lags, results.ens.mme_ratio, 'r.','MarkerSize',8);
end
yline(mr.bm_mme, 'k-', 'LineWidth',1.5);
if ~isnan(mr.ctrw_mme), yline(mr.ctrw_mme,'m--','LineWidth',1.5); end
if ~isnan(mr.fbm_mme),  yline(mr.fbm_mme, 'g--','LineWidth',1.5); end
if ~isnan(mr.obs_mme),  yline(mr.obs_mme, 'r-', 'LineWidth',1.5); end
% Bootstrap CI for obs_mme: horizontal shaded band
if hasBoot && ~any(isnan(boot.obs_mme_ci))
    xl3 = xlim; if xl3(1)==0 && xl3(2)==1, xl3=[lags(1) lags(end)]; end
    fill([xl3(1) xl3(2) xl3(2) xl3(1)], ...
         [boot.obs_mme_ci(1) boot.obs_mme_ci(1) boot.obs_mme_ci(2) boot.obs_mme_ci(2)], ...
         [1 0 0],'FaceAlpha',0.12,'EdgeColor','none');
end
xl = xlim;
text(xl(2)*0.55, mr.bm_mme+0.06,'BM (1.49)','Color','k','FontSize',8);
if ~isnan(mr.ctrw_mme)
    text(xl(2)*0.55,mr.ctrw_mme+0.06,sprintf('CTRW (%.2f)',mr.ctrw_mme),'Color','m','FontSize',8);
end
if ~isnan(mr.fbm_mme)
    text(xl(2)*0.55,mr.fbm_mme+0.06,sprintf('FBM (%.2f)',mr.fbm_mme),'Color','g','FontSize',8);
end
if ~isnan(mr.obs_mme)
    text(xl(2)*0.55,mr.obs_mme+0.06,sprintf('Obs (%.2f)',mr.obs_mme),'Color','r','FontSize',8);
end
ylim([0 max(4,(mr.ctrw_mme+0.5)*~isnan(mr.ctrw_mme)+4*isnan(mr.ctrw_mme))]);
xlabel('Lag (frames)'); ylabel('\langler^4_{max}\rangle/\langler^2_{max}\rangle^2');
title('MME Moment Ratio');

%-- Panel 4: Growing sphere -------------------------------------------------
subplot(2,3,4); hold on; box on; grid on;
aval = results.alpha_MSD;
if ~isnan(aval)
    ok = isfinite(gs.prob) & gs.prob > 0;
    if any(ok)
        xv = lags(ok).^(aval/2);
        plot(xv, gs.prob(ok),'k.','MarkerSize',8,'DisplayName','Observed');
        if ~isnan(gs.slope)
            b   = mean(log(gs.prob(ok)) - gs.slope*log(lags(ok)));
            xf  = linspace(min(xv),max(xv),200);
            tf  = xf.^(2/aval);
            plot(xf, exp(gs.slope*log(tf)+b),'r-','LineWidth',1.5,...
                 'DisplayName',sprintf('slope=%.3f, d_f\\approx%.2f',gs.slope,gs.df));
        end
    end
end
xlabel('t^{\alpha/2}'); ylabel('Pr(r \leq r_0 t^{\alpha/2})');
title('Growing Sphere'); legend('Location','best','FontSize',8);

%-- Panel 5: Classification text -------------------------------------------
subplot(2,3,5); axis off;
lines = { ...
    sprintf('\\bfProcess:\\rm  %s',      cls.process), ...
    sprintf('Confidence: %s',           cls.confidence), ...
    ' ', ...
    sprintf('\\alpha_{MSD}         = %.3f', results.alpha_MSD), ...
    sprintf('\\alpha''_{MME}         = %.3f', results.alpha_MME), ...
    sprintf('\\alpha_{corr(FBM)}  = %.3f', results.alpha_corrected), ...
    ' ', ...
    sprintf('Obs. reg. ratio = %.2f  (BM=%.2f)', mr.obs_reg, mr.bm_reg), ...
    sprintf('Obs. MME ratio  = %.2f  (BM=%.2f)', mr.obs_mme, mr.bm_mme), ...
    sprintf('d_f (GS)  \\approx  %.2f', gs.df), ...
    ' ', ...
    sprintf('Noise correction: %s', nc.mode), ...
};
if ~strcmp(nc.mode,'none')
    lines{end+1} = sprintf('  B_{MSD} = %.3g, B_{MME2} = %.3g', nc.B_MSD_ens, nc.B_MME2_ens);
    lines{end+1} = sprintf('  \\sigma_{loc} est. = %.3g', nc.sigma_loc);
end
% Bootstrap CIs
if hasBoot
    lines{end+1} = ' ';
    lines{end+1} = sprintf('Bootstrap 95%% CI  (n=%d)', boot.nBoot);
    lines{end+1} = sprintf('  \\alpha_{MSD}:  [%.3f, %.3f]', boot.alpha_MSD_ci(1), boot.alpha_MSD_ci(2));
    lines{end+1} = sprintf('  \\alpha''_{MME}: [%.3f, %.3f]', boot.alpha_MME_ci(1), boot.alpha_MME_ci(2));
    if ~any(isnan(boot.alpha_corr_ci))
        lines{end+1} = sprintf('  \\alpha_{corr}: [%.3f, %.3f]', boot.alpha_corr_ci(1), boot.alpha_corr_ci(2));
    end
    lines{end+1} = sprintf('  reg ratio:  [%.2f, %.2f]', boot.obs_reg_ci(1), boot.obs_reg_ci(2));
    lines{end+1} = sprintf('  MME ratio:  [%.2f, %.2f]', boot.obs_mme_ci(1), boot.obs_mme_ci(2));
    if ~any(isnan(boot.gs_slope_ci))
        lines{end+1} = sprintf('  GS slope:   [%.3f, %.3f]', boot.gs_slope_ci(1), boot.gs_slope_ci(2));
    end
    if ~any(isnan(boot.B_MSD_ci)) && ~strcmp(nc.mode,'none')
        lines{end+1} = sprintf('  B_{MSD} CI: [%.3g, %.3g]', boot.B_MSD_ci(1), boot.B_MSD_ci(2));
    end
end
for k = 1:min(numel(cls.notes),4)
    lines{end+1} = ['  \bullet ' cls.notes{k}]; %#ok<AGROW>
end
text(0.03, 0.97, lines, 'Units','normalized','VerticalAlignment','top', ...
     'FontSize',8,'Interpreter','tex');
title('Classification Summary');

%-- Panel 6: Effective sample size -----------------------------------------
subplot(2,3,6); hold on; box on; grid on;
hasEns = isfield(results,'ens');
hasTA  = isfield(results,'ta');
if hasEns && hasTA
    yyaxis left;
    plot(lags, results.ens.Neff,'b.-','DisplayName','N_{eff} (ens)');
    ylabel('N_{tracks} (ensemble)');
    yyaxis right;
    plot(lags, results.ta.Neff,'r.-','DisplayName','N_{win} (TA)');
    ylabel('N_{windows} (TA)');
elseif hasEns
    plot(lags, results.ens.Neff,'b.-');
    ylabel('N_{tracks} (ensemble)');
elseif hasTA
    plot(lags, results.ta.Neff,'r.-');
    ylabel('N_{windows} (TA)');
end
xlabel('Lag (frames)');
title('Effective sample size per lag');
legend('Location','best','FontSize',8);

end  % makePlots


%% =========================================================================
%%  HELPERS
%% =========================================================================

% addFitLine: draw K*(lag^alpha + B) on current log-log axes.
% The raw data has the offset still in it, so the plotted fit line is
% K*lag^alpha + B, which visually shows the offset correction.
function addFitLine(lags, vals_fit, alpha, B, linespec, name)
if isnan(alpha), return; end
ok = isfinite(vals_fit) & vals_fit > 0;
if sum(ok) < 2, return; end
% Estimate K from the corrected (offset-removed) data
K = exp(mean(log(vals_fit(ok)) - alpha * log(lags(ok))));
% Plot the fit curve shifted back by B so it tracks the raw data
loglog(lags, K * lags.^alpha + B, linespec, 'LineWidth',1.5, 'DisplayName',name);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

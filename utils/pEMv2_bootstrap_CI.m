function ciResults = pEMv2_bootstrap_CI(deltaX, results, trackInfo, params, varargin)
%PEMV2_BOOTSTRAP_CI  Non-parametric bootstrap confidence intervals for pEMv2.
%
%   ciResults = pEMv2_bootstrap_CI(deltaX, results, trackInfo, params)
%   ciResults = pEMv2_bootstrap_CI(..., 'nBoot', 500, 'alpha', 0.05, ...
%       'nRandomStarts', 5, 'covModel', 'fBM', 'verbose', true)
%
%   Inputs (all available in workspace after main_pEMv2.m runs):
%     deltaX    - cell array of displacement vectors
%     results   - pEMv2 output struct from pEMv2_SPT
%     trackInfo - struct with .vacf_exp, .numberOfTracks, .numFeatures,
%                 .splitLength, .dimensions, .lambda, .dt, .R
%     params    - struct with .converged, .maxiter
%
%   Name-value options:
%     'nBoot'         - number of bootstrap resamples (default 500)
%     'alpha'         - significance level for CI (default 0.05 = 95% CI)
%     'nRandomStarts' - random restarts per bootstrap replicate (default 5)
%     'covModel'      - covariance model: 'fBM' (default), 'normal', 'confined'
%     'verbose'       - display progress (default true)
%     'parallel'      - use parfor for bootstrap loop (default false,
%                       requires Parallel Computing Toolbox)
%     'resampleLevel' - 'subtracks' (default, matches pEM.m) or 'parent'
%                       (resample parent tracks, include all their
%                       sub-tracks; requires trackInfo.splitIndex).
%                       'parent' gives wider, more realistic CIs by
%                       respecting within-track correlation.
%
%   Output: ciResults struct with fields P, D, sigma, vacf (and alpha_fbm
%   for fBM or L for confined). Each contains .point, .ci_lo, .ci_hi,
%   .boot_dist.

% Parse name-value options
p = inputParser;
addParameter(p, 'nBoot', 500, @isnumeric);
addParameter(p, 'alpha', 0.05, @isnumeric);
addParameter(p, 'nRandomStarts', 5, @isnumeric);
addParameter(p, 'covModel', 'fBM', @ischar);
addParameter(p, 'verbose', true, @islogical);
addParameter(p, 'parallel', false, @islogical);
addParameter(p, 'resampleLevel', 'subtracks', ...
    @(x) ismember(x, {'subtracks','parent'}));
parse(p, varargin{:});
opts = p.Results;

% Extract point estimates
K = results.optimalSize;
vacf_opt = results.optimalVacf;   % K x numFeatures
P_opt = results.optimalP;         % 1 x K

nTracks = trackInfo.numberOfTracks;
numFeatures = trackInfo.numFeatures;
dt = trackInfo.dt;
R = trackInfo.R;

% Sort point estimates by vacf(:,1) for consistent ordering
[~, ord] = sort(vacf_opt(:,1));
vacf_opt = vacf_opt(ord,:);
P_opt = P_opt(ord);

% Compute point-estimate hard-assignment fractions
posteriorProb = results.posteriorProb;  % nTracks x K
posteriorProb = posteriorProb(:,ord);   % apply same sort order
[~, assignment] = max(posteriorProb, [], 2);
point_hardP = zeros(1, K);
for k = 1:K
    point_hardP(k) = sum(assignment == k) / nTracks;
end

% Fit point-estimate biophysical parameters
point_D = zeros(1, K);
point_sigma = zeros(1, K);
point_extra = zeros(1, K);  % alpha (fBM) or L (confined)
for k = 1:K
    C_exp = toeplitz(vacf_opt(k,:));
    fitP = make_fit_params(C_exp, dt, R);
    [pFit, ~] = OptimalParameters(C_exp, opts.covModel, fitP);
    point_D(k) = pFit(1);
    point_sigma(k) = abs(pFit(end));
    if length(pFit) == 3
        point_extra(k) = pFit(2);
    end
end

% Preallocate bootstrap storage
boot_P     = NaN(opts.nBoot, K);
boot_vacf  = NaN(opts.nBoot, K, numFeatures);
boot_D     = NaN(opts.nBoot, K);
boot_sigma = NaN(opts.nBoot, K);
boot_extra = NaN(opts.nBoot, K);  % alpha or L
boot_hardP = NaN(opts.nBoot, K);
valid      = false(opts.nBoot, 1);

% EM parameters (only converged and maxiter needed by EM.m)
emParams.converged = params.converged;
emParams.maxiter   = params.maxiter;

% Precompute parent-track groups for parent-level resampling
if strcmp(opts.resampleLevel, 'parent')
    if ~isfield(trackInfo, 'splitIndex')
        error('pEMv2_bootstrap_CI:noSplitIndex', ...
            'trackInfo.splitIndex is required for resampleLevel=''parent''. Set in main_pEMv2.m.');
    end
    parentGroups = accumarray(trackInfo.splitIndex, (1:nTracks)', [], @(x){x});
else
    parentGroups = {};  % unused, but needed as broadcast variable for parfor
end

% Launch parallel pool if requested
useParallel = opts.parallel;
if useParallel
    if isempty(ver('parallel'))
        warning('pEMv2_bootstrap_CI:noToolbox', ...
            'Parallel Computing Toolbox not found. Falling back to serial.');
        useParallel = false;
    elseif isempty(gcp('nocreate'))
        parpool;
    end
end

% Copy broadcast variables for parfor compatibility
covModel_b = opts.covModel;
nRandomStarts_b = opts.nRandomStarts;
resampleLevel_b = opts.resampleLevel;
verbose_b = opts.verbose && ~useParallel;  % suppress per-iteration output in parfor

% Bootstrap loop (parfor or for)
nBoot_ = opts.nBoot;
doVerbose = opts.verbose;
if useParallel
    % DataQueue for parfor progress tracking
    if doVerbose
        dq = parallel.pool.DataQueue;
        afterEach(dq, @(~) update_parfor_progress(nBoot_));
        fprintf('Bootstrap   0/%d', nBoot_);
    else
        dq = parallel.pool.DataQueue;  % dummy, never sent to
    end
    parfor b = 1:nBoot_
        [boot_P(b,:), boot_vacf(b,:,:), boot_D(b,:), boot_sigma(b,:), ...
            boot_extra(b,:), boot_hardP(b,:), valid(b)] = ...
            run_one_bootstrap(b, deltaX, trackInfo, vacf_opt, P_opt, ...
            emParams, K, nTracks, numFeatures, dt, R, ...
            nRandomStarts_b, covModel_b, resampleLevel_b, parentGroups, false);
        if doVerbose, send(dq, b); end
    end
    if doVerbose
        fprintf('\nBootstrap complete: %d/%d valid resamples\n', sum(valid), nBoot_);
    end
else
    for b = 1:opts.nBoot
        if verbose_b
            fprintf('Bootstrap %d/%d ...', b, opts.nBoot);
        end
        [boot_P(b,:), boot_vacf(b,:,:), boot_D(b,:), boot_sigma(b,:), ...
            boot_extra(b,:), boot_hardP(b,:), valid(b)] = ...
            run_one_bootstrap(b, deltaX, trackInfo, vacf_opt, P_opt, ...
            emParams, K, nTracks, numFeatures, dt, R, ...
            nRandomStarts_b, covModel_b, resampleLevel_b, parentGroups, verbose_b);
    end
end

% Extract valid samples
nGood     = sum(valid);
goodP     = boot_P(valid,:);
goodD     = boot_D(valid,:);
goodSigma = boot_sigma(valid,:);
goodExtra = boot_extra(valid,:);
goodHardP = boot_hardP(valid,:);
goodVacf  = boot_vacf(valid,:,:);

% Percentile CI bounds
lo = opts.alpha/2 * 100;
hi = (1 - opts.alpha/2) * 100;

% Assemble output struct
ciResults.P.point     = P_opt;
ciResults.P.ci_lo     = prctile(goodP, lo);
ciResults.P.ci_hi     = prctile(goodP, hi);
ciResults.P.boot_dist = goodP;

ciResults.hardP.point     = point_hardP;
ciResults.hardP.ci_lo     = prctile(goodHardP, lo);
ciResults.hardP.ci_hi     = prctile(goodHardP, hi);
ciResults.hardP.boot_dist = goodHardP;

ciResults.D.point     = point_D;
ciResults.D.ci_lo     = prctile(goodD, lo);
ciResults.D.ci_hi     = prctile(goodD, hi);
ciResults.D.boot_dist = goodD;

ciResults.sigma.point     = point_sigma;
ciResults.sigma.ci_lo     = prctile(goodSigma, lo);
ciResults.sigma.ci_hi     = prctile(goodSigma, hi);
ciResults.sigma.boot_dist = goodSigma;

ciResults.vacf.point     = vacf_opt;
ciResults.vacf.ci_lo     = squeeze(prctile(goodVacf, lo, 1));
ciResults.vacf.ci_hi     = squeeze(prctile(goodVacf, hi, 1));
ciResults.vacf.boot_dist = goodVacf;

if strcmp(opts.covModel, 'fBM')
    ciResults.alpha_fbm.point     = point_extra;
    ciResults.alpha_fbm.ci_lo     = prctile(goodExtra, lo);
    ciResults.alpha_fbm.ci_hi     = prctile(goodExtra, hi);
    ciResults.alpha_fbm.boot_dist = goodExtra;
end

if strcmp(opts.covModel, 'confined')
    ciResults.L.point     = point_extra;
    ciResults.L.ci_lo     = prctile(goodExtra, lo);
    ciResults.L.ci_hi     = prctile(goodExtra, hi);
    ciResults.L.boot_dist = goodExtra;
end

ciResults.covModel      = opts.covModel;
ciResults.resampleLevel = opts.resampleLevel;
ciResults.nBoot         = opts.nBoot;
ciResults.nGood         = nGood;
ciResults.alpha         = opts.alpha;

% Display results
if opts.verbose
    display_results(ciResults, K, opts);
end

% Warn if too many invalid resamples
if nGood < 0.9 * opts.nBoot
    warning('pEMv2_bootstrap_CI:lowValid', ...
        '%.1f%% of bootstrap resamples were invalid. States may not be well-separated.', ...
        (1 - nGood/opts.nBoot)*100);
end

end


%% Helper: run a single bootstrap replicate (parfor-compatible, no continue)
function [bP, bVacf, bD, bSigma, bExtra, bHardP, bValid] = ...
        run_one_bootstrap(~, deltaX, trackInfo, vacf_opt, P_opt, ...
        emParams, K, nTracks, numFeatures, dt, R, ...
        nRandomStarts, covModel, resampleLevel, parentGroups, verbose)

    % Default: invalid replicate (NaN storage)
    bP     = NaN(1, K);
    bVacf  = NaN(K, numFeatures);
    bD     = NaN(1, K);
    bSigma = NaN(1, K);
    bExtra = NaN(1, K);
    bHardP = NaN(1, K);
    bValid = false;

    % Resample
    if strcmp(resampleLevel, 'parent')
        % Resample parent tracks with replacement, include all sub-tracks
        nParents = length(parentGroups);
        pidx = ceil(rand(nParents,1) * nParents);
        idx = vertcat(parentGroups{pidx});
    else
        % Resample sub-tracks with replacement (matches pEM.m lines 25-31)
        idx = ceil(rand(nTracks,1) * nTracks);
    end
    nBoot = length(idx);
    bootDeltaX = deltaX(idx);
    bsInfo = trackInfo;
    bsInfo.vacf_exp = trackInfo.vacf_exp(idx,:,:);
    bsInfo.numberOfTracks = nBoot;

    % Warm start from point estimates
    try
        [v1, p1, g1, L1] = EM(bootDeltaX, vacf_opt, P_opt, emParams, bsInfo);
        bestV = v1; bestP = p1; bestG = g1; bestL = L1(end);
    catch
        if verbose, fprintf(' EM failed, skipping\n'); end
        return;
    end

    % Random restarts
    for r = 1:nRandomStarts
        try
            [v0, p0] = RandomInitialization(K, bsInfo.vacf_exp, 2);
            [vr, pr, gr, Lr] = EM(bootDeltaX, v0, p0, emParams, bsInfo);
            if Lr(end) > bestL
                bestV = vr; bestP = pr; bestG = gr; bestL = Lr(end);
            end
        catch
            % skip failed restart
        end
    end

    % Sort states by vacf(:,1) ascending to resolve label switching
    [~, sord] = sort(bestV(:,1));
    bestV = bestV(sord,:);
    bestP = bestP(sord);
    bestG = bestG(:,sord,:);

    % Degeneracy check
    if any(bestP < 1e-6) || any(isnan(bestV(:)))
        if verbose, fprintf(' degenerate, skipping\n'); end
        return;
    end

    % Fit covariance model for each state
    tmpD = zeros(1, K);
    tmpSigma = zeros(1, K);
    tmpExtra = zeros(1, K);
    for k = 1:K
        C_exp = toeplitz(bestV(k,:));
        fitP = make_fit_params(C_exp, dt, R);
        try
            [pFit, ~] = OptimalParameters(C_exp, covModel, fitP);
            tmpD(k) = pFit(1);
            tmpSigma(k) = abs(pFit(end));
            if length(pFit) == 3
                tmpExtra(k) = pFit(2);
            end
        catch
            if verbose, fprintf(' fit failed, skipping\n'); end
            return;
        end
    end

    % Compute hard-assignment fractions from posterior
    dim = size(bestG, 3);
    postProb = sum(bestG, 3) / dim;  % nBoot x K (matches pEM.m line 64)
    [~, asgn] = max(postProb, [], 2);
    tmpHardP = zeros(1, K);
    for k = 1:K
        tmpHardP(k) = sum(asgn == k) / nBoot;
    end

    % All steps succeeded — store results
    bP     = bestP;
    bVacf  = bestV;
    bD     = tmpD;
    bSigma = tmpSigma;
    bExtra = tmpExtra;
    bHardP = tmpHardP;
    bValid = true;

    if verbose, fprintf(' done\n'); end
end


%% Helper: parfor progress counter (called via DataQueue)
function update_parfor_progress(nBoot)
    persistent count
    if isempty(count), count = 0; end
    count = count + 1;
    fprintf(repmat('\b', 1, length(sprintf('%d/%d', count-1, nBoot))));
    fprintf('%d/%d', count, nBoot);
    if count == nBoot, count = []; end  % reset for next call
end


%% Helper: construct initial parameters for OptimalParameters
function fitP = make_fit_params(C_exp, dt, R)
    fitP.dt     = dt;
    fitP.R      = R;
    fitP.D0     = (C_exp(1,1) + 2*C_exp(1,2)) / (2*dt);
    fitP.sigma0 = sqrt(abs(C_exp(1,1)/2 - fitP.D0*dt*(1-2*R)));
    fitP.L0     = fitP.D0 * 5;
    fitP.A0     = 1;
end


%% Helper: display formatted CI table
function display_results(ci, K, opts)
    fprintf('\n==================================================================\n');
    fprintf('BOOTSTRAP CONFIDENCE INTERVALS (%d%% CI, %d/%d valid resamples)\n', ...
        round((1-opts.alpha)*100), ci.nGood, opts.nBoot);
    fprintf('Covariance model: %s  |  Resample level: %s\n', ci.covModel, ci.resampleLevel);
    fprintf('==================================================================\n');
    for k = 1:K
        fprintf('\nState %d:\n', k);
        fprintf('  P     = %.4f  [%.4f, %.4f]  (soft)\n', ...
            ci.P.point(k), ci.P.ci_lo(k), ci.P.ci_hi(k));
        fprintf('  P     = %.4f  [%.4f, %.4f]  (hard assignment)\n', ...
            ci.hardP.point(k), ci.hardP.ci_lo(k), ci.hardP.ci_hi(k));
        fprintf('  D     = %.4f  [%.4f, %.4f] um^2/s\n', ...
            ci.D.point(k), ci.D.ci_lo(k), ci.D.ci_hi(k));
        if strcmp(opts.covModel, 'fBM')
            fprintf('  alpha = %.4f  [%.4f, %.4f]\n', ...
                ci.alpha_fbm.point(k), ci.alpha_fbm.ci_lo(k), ci.alpha_fbm.ci_hi(k));
        end
        if strcmp(opts.covModel, 'confined')
            fprintf('  L     = %.4f  [%.4f, %.4f] um\n', ...
                ci.L.point(k), ci.L.ci_lo(k), ci.L.ci_hi(k));
        end
        fprintf('  sigma = %.4f  [%.4f, %.4f] um\n', ...
            ci.sigma.point(k), ci.sigma.ci_lo(k), ci.sigma.ci_hi(k));
    end
    fprintf('\n==================================================================\n');
end

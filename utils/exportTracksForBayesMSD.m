function exportTracksForBayesMSD(tc, out_path, varargin)
%EXPORTTRACKSFORBAYESMSD  Save tracks from a TrajectoryCollection for bayesmsd.
%
%   exportTracksForBayesMSD(tc, '/path/to/out.mat')
%   exportTracksForBayesMSD(tc, '/path/to/out.mat', 'TrackType', 'culled', ...)
%
%   Saves a .mat file containing:
%     tracks    - cell array (N x 1) of (L x 2) double matrices [x µm, y µm]
%     dt        - frame interval in seconds
%     frac      - ExposureFraction (Te/dt) used for fBM fitting
%     n_tracks  - number of tracks saved
%     track_type - string: 'culled' | 'raw' | 'relative'
%     matlab_K     - population MLE K  (if FBMResults available)
%     matlab_alpha - population MLE alpha
%     matlab_sigma - population MLE sigma
%
%   Options:
%     'TrackType'        - 'culled' (default) | 'raw' | 'relative'
%     'Condition'        - restrict to one condition (default: all)
%     'MinTrackLength'   - minimum frames (default: 15)
%     'MaxTracks'        - maximum tracks to save; randomly subsample if
%                          more are available (default: Inf = save all)
%     'MinStepVar'       - minimum mean squared step size µm² (default: 0)
%     'ExposureFraction' - Te/dt stored as metadata (default: 1)
%     'RLState'          - if set, save only tracks from this RL state index.
%                          Requires tc.getRLDecomposition() to have been run.
%     'Seed'             - random seed for subsampling (default: 42)

p = inputParser;
addRequired(p,  'tc',       @(x) isa(x, 'TrajectoryCollection'));
addRequired(p,  'out_path', @ischar);
addParameter(p, 'TrackType',        'culled', @(x) ismember(lower(char(x)), {'culled','raw','relative'}));
addParameter(p, 'Condition',        [],       @(x) ischar(x)||isstring(x)||iscell(x)||isempty(x));
addParameter(p, 'MinTrackLength',   15,       @isnumeric);
addParameter(p, 'MaxTracks',        Inf,      @isnumeric);
addParameter(p, 'MinStepVar',       0,        @isnumeric);
addParameter(p, 'ExposureFraction', 1,        @isnumeric);
addParameter(p, 'dt',               NaN,      @isnumeric);
addParameter(p, 'RLState',          [],       @isnumeric);
addParameter(p, 'Seed',             42,       @isnumeric);
parse(p, tc, out_path, varargin{:});
o = p.Results;

trackType = lower(char(o.TrackType));

% ── Retrieve tracks ──────────────────────────────────────────────────────────

if ~isempty(o.RLState)
    % Pull a specific RL-classified state
    if ~isfield(tc.RLResults, 'computed_quantities') || isempty(tc.RLResults.computed_quantities)
        error('exportTracksForBayesMSD: run tc.getRLDecomposition() first to use RLState.');
    end
    cq      = tc.RLResults.computed_quantities;
    idx_map = tc.RLResults.track_index_map;
    rl_type = tc.RLResults.track_type;
    s       = round(o.RLState);
    if s < 1 || s > numel(cq.classified_tracks)
        error('exportTracksForBayesMSD: RLState %d out of range (1–%d).', s, numel(cq.classified_tracks));
    end
    rl_idx = cq.classified_tracks{s};
    idx    = idx_map(rl_idx);

    switch rl_type
        case 'relative', all_tracks = tc.getAllRelativeTracks();
        case 'raw',      all_tracks = tc.getAllRawTracks();
        otherwise,       all_tracks = tc.getAllCulledTracks();
    end
    tracks = all_tracks(idx);
    trackType = rl_type;
    fprintf('exportTracksForBayesMSD: RL state %d — %d tracks\n', s, numel(tracks));
else
    switch trackType
        case 'culled'
            if isempty(o.Condition)
                tracks = tc.getAllCulledTracks();
            else
                tracks = tc.getTracksByCondition(o.Condition);
            end
        case 'relative'
            if isempty(o.Condition)
                tracks = tc.getAllRelativeTracks();
            else
                tracks = tc.getRelativeTracksByCondition(o.Condition);
            end
        case 'raw'
            if isempty(o.Condition)
                tracks = tc.getAllRawTracks();
            else
                tracks = tc.getTracksByCondition(o.Condition, 'UseRaw', true);
            end
    end
end

if isempty(tracks)
    error('exportTracksForBayesMSD: no tracks found (TrackType=''%s'').', trackType);
end

% ── Filter ───────────────────────────────────────────────────────────────────

Lmin = round(o.MinTrackLength);
keep = cellfun(@(t) size(t,1) >= Lmin, tracks);
tracks = tracks(keep);
fprintf('exportTracksForBayesMSD: %d tracks after length filter (>= %d frames)\n', numel(tracks), Lmin);

if o.MinStepVar > 0
    mobile = cellfun(@(t) mean(diff(t(:,1)).^2 + diff(t(:,2)).^2) >= o.MinStepVar, tracks);
    tracks = tracks(mobile);
    fprintf('exportTracksForBayesMSD: %d tracks after MinStepVar filter\n', numel(tracks));
end

if isempty(tracks)
    error('exportTracksForBayesMSD: no tracks remain after filtering.');
end

% ── Subsample ────────────────────────────────────────────────────────────────

N = numel(tracks);
if isfinite(o.MaxTracks) && N > o.MaxTracks
    rng(o.Seed);
    idx   = randperm(N, round(o.MaxTracks));
    tracks = tracks(sort(idx));
    fprintf('exportTracksForBayesMSD: subsampled to %d tracks (seed=%d)\n', numel(tracks), o.Seed);
end

% ── Strip to [x, y] only and ensure double ───────────────────────────────────

for k = 1:numel(tracks)
    tracks{k} = double(tracks{k}(:, 1:2));
end

% ── Metadata ─────────────────────────────────────────────────────────────────

% Resolve dt from stored results or explicit argument (private method not accessible)
dt   = NaN;
frac = o.ExposureFraction;

if ~isnan(o.dt)
    dt = o.dt;
elseif isstruct(tc.RLResults) && isfield(tc.RLResults, 'dt') && ~isempty(tc.RLResults.dt)
    dt   = tc.RLResults.dt;
    frac = tc.RLResults.frac;
elseif isstruct(tc.FBMResults) && isfield(tc.FBMResults, 'dt') && ~isempty(tc.FBMResults.dt)
    dt = tc.FBMResults.dt;
elseif isstruct(tc.MSDResults) && isfield(tc.MSDResults, 'dt') && ~isempty(tc.MSDResults.dt)
    dt = tc.MSDResults.dt;
end

if isnan(dt)
    error(['exportTracksForBayesMSD: cannot determine dt automatically. ' ...
           'Pass it explicitly: exportTracksForBayesMSD(tc, path, ''dt'', 0.02)']);
end

matlab_K     = NaN;
matlab_alpha = NaN;
matlab_sigma = NaN;

if isstruct(tc.FBMResults) && isfield(tc.FBMResults, 'K')
    matlab_K     = tc.FBMResults.K;
    matlab_alpha = tc.FBMResults.alpha;
    matlab_sigma = tc.FBMResults.sigma;
end

n_tracks = numel(tracks);

% ── Save ─────────────────────────────────────────────────────────────────────

% Use -v7 (MATLAB default) for files < 2 GB; scipy.io.loadmat reads this fine.
% Use -v7.3 (HDF5) for larger files — requires h5py in Python.
file_size_est_mb = n_tracks * 50 * 2 * 8 / 1e6;  % rough estimate

if file_size_est_mb > 1800
    fmt = '-v7.3';
    fprintf('exportTracksForBayesMSD: large file (~%.0f MB est.), saving as HDF5 (v7.3)\n', file_size_est_mb);
    fprintf('  Python: requires h5py  (pip install h5py)\n');
else
    fmt = '-v7';
end

track_type = trackType; %#ok<NASGU>
save(out_path, 'tracks', 'dt', 'frac', 'n_tracks', 'track_type', ...
    'matlab_K', 'matlab_alpha', 'matlab_sigma', fmt);

fprintf('exportTracksForBayesMSD: saved %d tracks to %s\n', n_tracks, out_path);
fprintf('  dt=%.4f s, ExposureFraction=%.3f\n', dt, frac);
if ~isnan(matlab_K)
    fprintf('  MATLAB MLE: alpha=%.3f, K=%.4g µm²/s^alpha, sigma=%.4f µm\n', ...
        matlab_alpha, matlab_K, matlab_sigma);
end
fprintf('\nRun in Python:\n');
fprintf('  python run_bayesmsd.py "%s" \\\n', out_path);
fprintf('      --dt %.4f \\\n', dt);
fprintf('      --exposure_fraction %.3f', frac);
if ~isnan(matlab_alpha)
    fprintf(' \\\n      --matlab_K %.4g --matlab_alpha %.4f --matlab_sigma %.4f', ...
        matlab_K, matlab_alpha, matlab_sigma);
end
fprintf('\n');

end

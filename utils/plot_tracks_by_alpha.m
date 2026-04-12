function plot_tracks_by_alpha(tc, alpha_dist, varargin)
%PLOT_TRACKS_BY_ALPHA  Plot individual tracks coloured by their per-track alpha.
%
%   plot_tracks_by_alpha(tc, alpha_dist)
%   plot_tracks_by_alpha(tc, alpha_dist, 'AlphaMin', 1.0, 'N', 20)
%
%   Selects tracks whose fitted alpha falls in [AlphaMin, AlphaMax],
%   plots up to N of them, and colour-codes each by its alpha value.
%   Useful for visually checking whether high-alpha tracks are genuinely
%   superdiffusive or artefactual (many gaps, short segments, etc.).
%
%   Inputs:
%     tc          - TrajectoryCollection (must have been the source for alpha_dist)
%     alpha_dist  - struct returned by tc.getFBMAlphaDistribution()
%
%   Name-value options:
%     'AlphaMin'   - lower alpha cutoff (default: 1.0)
%     'AlphaMax'   - upper alpha cutoff (default: Inf)
%     'TrackType'  - 'culled' (default) | 'relative' | 'raw'
%     'N'          - max tracks to plot (default: 20)
%     'Seed'       - rng seed for random selection (default: 1)
%     'ShowGaps'   - mark gap positions with red dots (default: true)
%     'Title'      - figure title string (default: auto)

p = inputParser;
addRequired(p,  'tc');
addRequired(p,  'alpha_dist');
addParameter(p, 'AlphaMin',  1.0,      @isnumeric);
addParameter(p, 'AlphaMax',  Inf,      @isnumeric);
addParameter(p, 'TrackType', 'culled', @ischar);
addParameter(p, 'N',         20,       @isnumeric);
addParameter(p, 'Seed',      1,        @isnumeric);
addParameter(p, 'ShowGaps',  true,     @islogical);
addParameter(p, 'Title',     '',       @(x) ischar(x)||isstring(x));
parse(p, tc, alpha_dist, varargin{:});
o = p.Results;

% ── Retrieve the same track population used for alpha_dist ───────────
switch lower(o.TrackType)
    case 'culled'
        tracks = tc.getAllCulledTracks();
    case 'relative'
        tracks = tc.getAllRelativeTracks();
    case 'raw'
        tracks = tc.getAllRawTracks();
    otherwise
        error('TrackType must be culled, relative, or raw.');
end

% alpha_dist.original_index maps each result entry back to the
% corresponding position in the full culled/relative/raw track array.
% This is required because fitFBM_pertracks filters tracks by MinTrackLength
% and MinStepVar before fitting, so alpha_dist index i ≠ tracks{i}.
if ~isfield(alpha_dist, 'original_index')
    error(['alpha_dist does not contain original_index. ' ...
           'Re-run getFBMAlphaDistribution with ForceRecompute=true ' ...
           'to regenerate results with the updated fitFBM_pertracks.']);
end

% ── Find tracks in the selected alpha range ───────────────────────────
sel = find(alpha_dist.alpha >= o.AlphaMin & alpha_dist.alpha <= o.AlphaMax);
if isempty(sel)
    fprintf('No tracks found with alpha in [%.2f, %.2f].\n', o.AlphaMin, o.AlphaMax);
    return;
end

fprintf('Found %d tracks with alpha in [%.2f, %.2f].\n', ...
    numel(sel), o.AlphaMin, o.AlphaMax);

% Random subsample
rng(o.Seed);
if numel(sel) > o.N
    sel = sel(randperm(numel(sel), o.N));
end

% ── Colourmap: alpha → colour ─────────────────────────────────────────
all_alpha = alpha_dist.alpha(sel);
clim      = [min(all_alpha), max(all_alpha)];
if clim(1) == clim(2), clim(2) = clim(1) + 0.01; end
cmap      = jet(256);
cmap_idx  = @(a) max(1, min(256, round(255*(a - clim(1))/(clim(2) - clim(1))) + 1));

% ── Plot ──────────────────────────────────────────────────────────────
ttl = char(o.Title);
if isempty(ttl)
    ttl = sprintf('Tracks with \\alpha \\in [%.2f, %.2f]  (n=%d shown)', ...
        o.AlphaMin, o.AlphaMax, numel(sel));
end

figure('Name', ttl);
hold on;

for j = 1:numel(sel)
    idx      = sel(j);
    orig_idx = alpha_dist.original_index(idx);   % position in full track array
    a        = alpha_dist.alpha(idx);
    clr      = cmap(cmap_idx(a), :);

    if orig_idx <= numel(tracks)
        tr = tracks{orig_idx};
        x  = tr(:, 1);
        y  = tr(:, 2);

        plot(x, y, '-', 'Color', [clr 0.6], 'LineWidth', 1.2);
        plot(x(1), y(1), 'o', 'Color', clr, 'MarkerSize', 5, 'MarkerFaceColor', clr);

        % Show gap positions (frames where diff(frame) > 1)
        if o.ShowGaps && size(tr, 2) >= 3
            f    = tr(:, 3);
            gpos = find(diff(f) > 1);   % row index before each gap
            if ~isempty(gpos)
                plot(x(gpos), y(gpos), 'rx', 'MarkerSize', 8, 'LineWidth', 1.5);
            end
        end

        % Annotate with alpha value
        plot(x(end), y(end), 's', 'Color', clr, 'MarkerSize', 5, 'MarkerFaceColor', 'w');
        text(x(end), y(end), sprintf('  \\alpha=%.2f', a), ...
            'FontSize', 7, 'Color', clr * 0.7);
    end
end

% Colourbar
colormap(cmap);
cb = colorbar;
caxis(clim);
cb.Label.String = '\alpha';

axis equal; box off;
set(gca, 'FontSize', 13, 'FontWeight', 'bold');
xlabel('x  (\mum)');
ylabel('y  (\mum)');
title(ttl);

if o.ShowGaps
    % Legend entry for gap markers
    plot(NaN, NaN, 'rx', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Gap');
    legend('show', 'Location', 'best', 'Box', 'off');
end

hold off;

fprintf('Plotted %d tracks. Red x = gap position (frame discontinuity).\n', numel(sel));
end

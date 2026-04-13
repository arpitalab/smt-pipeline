%% ANALYSIS_STIFFNESS_COMPARISON
%  Composite analysis comparing H2B and ER dynamics on 1 kPa vs 100 kPa substrates.
%
%  Requires in workspace:
%    H2B_1k_E2, H2B_100k_E2   — TrajectoryCollection objects
%    ER_1k_slow, ER_100k_slow  — TrajectoryCollection objects
%  All must have fBM (getFBMParameters) and RL (getRLDecomposition) already computed.
%
%  Produces:
%    Figure 1 — MSD + fBM model curves (2 x 2)
%    Figure 2 — fBM alpha and Ka per RL state per condition
%    Figure 3 — Packing fraction histograms per RL state
%    Figure 4 — DACF with noise-corrected theory

% ── Parameters ────────────────────────────────────────────────────────────────
dt          = 0.2;          % s (frame interval)
sigma_fixed = 0.040;        % µm (localization precision)
thresh      = 4*sigma_fixed^2;
frac        = 0.05;         % ExposureFraction (Te/dt)
W_pc        = 80;           % packing fraction window (frames)
maxlag_dacf = 10;           % lags for DACF

% ── Collections ───────────────────────────────────────────────────────────────
tcs    = {H2B_1k_E2,     H2B_100k_E2,     ER_1k_slow,      ER_100k_slow};
labels = {'H2B 1 kPa',   'H2B 100 kPa',   'ER 1 kPa',      'ER 100 kPa'};
clrs   = {[0.25 0.50 0.85], [0.05 0.20 0.60], [0.90 0.45 0.10], [0.65 0.15 0.00]};
n_cond = numel(tcs);

% ── 1. Per-RL-state fBM (run once; uses cache if already computed) ─────────────
fprintf('=== getFBMByRLState ===\n');
fbm_states = cell(n_cond, 1);
for c = 1:n_cond
    fbm_states{c} = tcs{c}.getFBMByRLState( ...
        'CIMethod',      'profile',  ...
        'SubtrackLength', 20,        ...
        'MinStepVar',     thresh,    ...
        'MinAlpha',       0.1,       ...
        'PerTrack',       true,      ...
        'Verbose',        false);
end

% ── 2. Packing fraction per RL state ─────────────────────────────────────────
fprintf('=== computePackingFraction ===\n');
pf = cell(n_cond, 1);
for c = 1:n_cond
    cq      = tcs{c}.RLResults.computed_quantities;
    idx_map = tcs{c}.RLResults.track_index_map;
    tracks  = tcs{c}.getAllCulledTracks();
    n_states = numel(cq.classified_tracks);
    pf{c}    = cell(n_states, 1);
    for s = 1:n_states
        idx          = idx_map(cq.classified_tracks{s});
        state_tracks = tracks(idx);
        pf{c}{s}     = computePackingFraction(state_tracks, ...
            'WindowSize',    W_pc,  ...
            'MinTrackLength', W_pc+1, ...
            'dt',            dt);
        fprintf('  %s  state %d: %d tracks with pc\n', labels{c}, s, numel(pf{c}{s}.pc));
    end
end

% ── 3. DACF ───────────────────────────────────────────────────────────────────
fprintf('=== computeDACF ===\n');
dacf_vals = cell(n_cond, 1);
for c = 1:n_cond
    [dacf_vals{c}, lags_dacf] = tcs{c}.computeDACF( ...
        'TrackType',     'relative', ...
        'MaxLag',         maxlag_dacf, ...
        'MinTrackLength', 20,          ...
        'MinStepVar',     thresh);
end

% =============================================================================
%  FIGURE 1 — MSD + fBM model curves
% =============================================================================
figure('Name', 'MSD + fBM fit', 'Position', [50 50 1100 800]);
protein_pairs = {[1 2], [3 4]};   % {H2B pair, ER pair}
protein_names = {'H2B', 'ER'};

for pp = 1:2
    c1 = protein_pairs{pp}(1);
    c2 = protein_pairs{pp}(2);
    ax = subplot(1, 2, pp);
    hold on;

    for c = [c1, c2]
        if ~isfield(tcs{c}.MSDResults,'ensemble_mean') || isempty(tcs{c}.MSDResults.ensemble_mean), continue; end
        msd  = tcs{c}.MSDResults;
        % lag_axis may be absent in older cached results — reconstruct from lag_frames
        if isfield(msd, 'lag_axis') && ~isempty(msd.lag_axis)
            lags = msd.lag_axis;
        else
            lags = msd.lag_frames * msd.dt;
        end
        y    = msd.ensemble_mean;

        % Empirical MSD
        plot(lags, y, '-', 'Color', clrs{c}, 'LineWidth', 2, 'DisplayName', labels{c});

        % fBM model curve (boot_fits median: [G, sigma_sq, alpha])
        if isfield(msd, 'boot_fits') && ~isempty(msd.boot_fits)
            med_fp   = median(msd.boot_fits, 1);
            lag_fr   = msd.lag_frames;
            frac_msd = msd.frac;
            dt_msd   = msd.dt;
            ycurve   = local_msd_model(med_fp, lag_fr, dt_msd, frac_msd);
            plot(lags, ycurve, '--', 'Color', clrs{c}*0.6 + 0.4, 'LineWidth', 1.2, ...
                'HandleVisibility', 'off');
        end
    end

    set(ax, 'XScale', 'log', 'YScale', 'log', 'FontSize', 12);
    xlabel('Lag time (s)'); ylabel('MSD (µm²)');
    title(protein_names{pp});
    legend('Location', 'northwest', 'Box', 'off');
    box off;
end
sgtitle('MSD: dashed = fBM model (median bootstrap fit)', 'FontSize', 13);

% =============================================================================
%  FIGURE 2 — fBM alpha and Ka per RL state
% =============================================================================
figure('Name', 'fBM parameters per RL state', 'Position', [50 100 1200 500]);

% Collect alpha and Ka for each condition x state
max_states = max(cellfun(@numel, fbm_states));

ax1 = subplot(1, 2, 1);  hold on;
ax2 = subplot(1, 2, 2);  hold on;

x_pos = 0;
xtick_pos = [];
xtick_lbl = {};

for c = 1:n_cond
    for s = 1:2
        st = fbm_states{c}(s);
        if st.skipped, continue; end
        x_pos = x_pos + 1;

        % Alpha with CI
        ci_a = st.CI.alpha;
        errorbar(ax1, x_pos, st.alpha, st.alpha - ci_a(1), ci_a(2) - st.alpha, ...
            'o', 'Color', clrs{c}, 'MarkerFaceColor', clrs{c}, ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 6);

        % Ka with CI
        ci_k = st.CI.Ka;
        errorbar(ax2, x_pos, st.Ka, st.Ka - ci_k(1), ci_k(2) - st.Ka, ...
            'o', 'Color', clrs{c}, 'MarkerFaceColor', clrs{c}, ...
            'MarkerSize', 8, 'LineWidth', 1.5, 'CapSize', 6);

        xtick_pos(end+1) = x_pos; %#ok<SAGROW>
        xtick_lbl{end+1} = sprintf('%s\nS%d', strrep(labels{c},' ','\n'), s); %#ok<SAGROW>
    end
    x_pos = x_pos + 0.5;   % gap between conditions
end

% Alpha panel
yline(ax1, 0.5,  '--k', 'Alpha', 0.4, 'LineWidth', 1);
yline(ax1, 0.40, ':k',  'Alpha', 0.4, 'LineWidth', 1);
set(ax1, 'XTick', xtick_pos, 'XTickLabel', xtick_lbl, 'XTickLabelRotation', 0, ...
    'FontSize', 11, 'TickLabelInterpreter', 'tex');
ylabel(ax1, '\alpha');
title(ax1, 'Anomalous exponent per state');
text(ax1, ax1.XLim(2), 0.5,  ' Rouse',          'FontSize', 9, 'Color', [0.4 0.4 0.4]);
text(ax1, ax1.XLim(2), 0.40, ' Fractal globule', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
box(ax1, 'off');

% Ka panel
set(ax2, 'XTick', xtick_pos, 'XTickLabel', xtick_lbl, 'XTickLabelRotation', 0, ...
    'YScale', 'log', 'FontSize', 11, 'TickLabelInterpreter', 'tex');
ylabel(ax2, 'K\alpha  (µm² s^{-\alpha})');
title(ax2, 'Transport coefficient per state');
box(ax2, 'off');

sgtitle('fBM parameters: error bars = 95% profile CI', 'FontSize', 13);

% Dummy legend
for c = 1:n_cond
    plot(ax1, NaN, NaN, 'o', 'Color', clrs{c}, 'MarkerFaceColor', clrs{c}, ...
        'MarkerSize', 8, 'DisplayName', labels{c});
end
legend(ax1, 'Location', 'best', 'Box', 'off');

edges = linspace(-2, 4, 40);   % log10(pc) edges
xf    = linspace(min(edges), max(edges), 200);

%% =============================================================================
%  FIGURE 3a — pc histograms: stiffness comparison within protein
%              rows = states 1 & 2,  cols = H2B | ER
% =============================================================================
figure('Name', 'PC: stiffness comparison', 'Position', [100 50 1100 700]);

for s = 1:2
    % H2B: 1kPa vs 100kPa
    subplot(2, 2, (s-1)*2 + 1);  hold on;
    for c = 1:2
        if s > numel(pf{c}), continue; end
        v = log10(pf{c}{s}.pc);  v = v(isfinite(v));
        histogram(v, edges, 'Normalization', 'pdf', ...
            'FaceColor', clrs{c}, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
            'DisplayName', labels{c});
        [mu, sg] = normfit(v);
        plot(xf, normpdf(xf, mu, sg), '-', 'Color', clrs{c}, 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end
    xlabel('log_{10}(packing fraction)'); ylabel('PDF');
    title(sprintf('H2B — state %d  (1 kPa vs 100 kPa)', s));
    legend('Box', 'off', 'Location', 'best');  box off;

    % ER: 1kPa vs 100kPa
    subplot(2, 2, (s-1)*2 + 2);  hold on;
    for c = 3:4
        if s > numel(pf{c}), continue; end
        v = log10(pf{c}{s}.pc);  v = v(isfinite(v));
        histogram(v, edges, 'Normalization', 'pdf', ...
            'FaceColor', clrs{c}, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
            'DisplayName', labels{c});
        [mu, sg] = normfit(v);
        plot(xf, normpdf(xf, mu, sg), '-', 'Color', clrs{c}, 'LineWidth', 1.8, ...
            'HandleVisibility', 'off');
    end
    xlabel('log_{10}(packing fraction)'); ylabel('PDF');
    title(sprintf('ER — state %d  (1 kPa vs 100 kPa)', s));
    legend('Box', 'off', 'Location', 'best');  box off;
end
sgtitle('Packing fraction: stiffness effect within protein', 'FontSize', 13);

%% =============================================================================
%  FIGURE 3b — pc histograms: H2B vs ER at same stiffness
%              rows = states 1 & 2,  cols = 1 kPa | 100 kPa
% =============================================================================
% Map: condition index → stiffness label and H2B/ER pairing
% H2B_1k=1, H2B_100k=2, ER_1k=3, ER_100k=4
stiff_pairs  = {[1 3], [2 4]};          % {1kPa pair, 100kPa pair}
stiff_labels = {'1 kPa', '100 kPa'};

figure('Name', 'PC: H2B vs ER', 'Position', [150 100 1100 700]);

for s = 1:2
    for sp = 1:2
        c_h2b = stiff_pairs{sp}(1);
        c_er  = stiff_pairs{sp}(2);

        subplot(2, 2, (s-1)*2 + sp);  hold on;

        pair_idx  = [c_h2b, c_er];
        pair_lbls = {'H2B', 'ER'};
        for pi = 1:2
            c   = pair_idx(pi);
            lbl = pair_lbls{pi};
            if s > numel(pf{c}), continue; end
            v = log10(pf{c}{s}.pc);  v = v(isfinite(v));
            histogram(v, edges, 'Normalization', 'pdf', ...
                'FaceColor', clrs{c}, 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
                'DisplayName', lbl);
            [mu, sg] = normfit(v);
            plot(xf, normpdf(xf, mu, sg), '-', 'Color', clrs{c}, 'LineWidth', 1.8, ...
                'HandleVisibility', 'off');
        end
        xlabel('log_{10}(packing fraction)'); ylabel('PDF');
        title(sprintf('%s — state %d', stiff_labels{sp}, s));
        legend('Box', 'off', 'Location', 'best');  box off;
    end
end
sgtitle('Packing fraction: H2B vs ER at same stiffness', 'FontSize', 13);

%% =============================================================================
%  FIGURE 4 — DACF with noise-corrected theory
% =============================================================================
figure('Name', 'DACF', 'Position', [150 150 600 420]);
hold on;

for c = 1:n_cond
    plot(lags_dacf, dacf_vals{c}, '-o', 'Color', clrs{c}, 'LineWidth', 2, ...
        'MarkerSize', 5, 'DisplayName', labels{c});
end

% Noise-corrected theory for each condition
for c = 1:n_cond
    fbm = tcs{c}.FBMResults;
    if isempty(fieldnames(fbm)), continue; end
    nc   = local_dacf_theory(fbm, maxlag_dacf);
    plot(lags_dacf, nc, '--', 'Color', clrs{c}, 'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end

yline(0, 'Color', [0.7 0.7 0.7], ':');
xlabel('Lag (frames)'); ylabel('DACF');
set(gca, 'FontSize', 12);
legend('Box', 'off', 'Location', 'best');
title('DACF: dashed = noise-corrected fBM theory');
box off;

% =============================================================================
%  SUMMARY TABLE
% =============================================================================
fprintf('\n%s\n', repmat('=',1,80));
fprintf('%-20s  %-6s  %-8s  %-10s  %-8s  %-8s\n', ...
    'Condition', 'State', 'N tracks', 'alpha', 'Ka', 'sigma');
fprintf('%s\n', repmat('-',1,80));
for c = 1:n_cond
    for s = 1:numel(fbm_states{c})
        st = fbm_states{c}(s);
        if st.skipped, continue; end
        fprintf('%-20s  %-6d  %-8d  %-10.3f  %-8.4g  %-8.4g\n', ...
            labels{c}, s, st.n_tracks, st.alpha, st.Ka, st.sigma);
    end
end
fprintf('%s\n', repmat('=',1,80));

% =============================================================================
%  Local functions
% =============================================================================

function curve = local_msd_model(fp, lag_frames, dt, frac)
% Reconstruct fBM MSD model from boot_fits parameters [G, sigma_sq, alpha].
G = fp(1); sig2 = fp(2); alpha = fp(3);
r     = frac ./ lag_frames;
b     = (abs(1 + r).^(2+alpha) + abs(1-r).^(2+alpha) - 2) ./ r.^2;
curve = G / ((1+alpha)*(2+alpha)) .* ((dt.*lag_frames).^alpha .* b - ...
        2*(frac*dt)^alpha) + 2*sig2;
end


function dacf_nc = local_dacf_theory(fbm, maxlag)
% Noise-corrected theoretical DACF from getFBMParameters output struct.
if ~isfield(fbm, 'alpha') || isempty(fbm.alpha)
    dacf_nc = nan(maxlag, 1);
    return;
end
K_m   = fbm.K;
a_m   = fbm.alpha;
s_m   = fbm.sigma;
dt_m  = fbm.dt;
Te    = fbm.frac * dt_m;

psi = @(tau) ((tau + Te).^(a_m+2) + abs(tau - Te).^(a_m+2) - 2*tau.^(a_m+2)) ...
             ./ (Te^2 * (a_m+1) * (a_m+2));

cov_lag = @(m) K_m * (psi((m+1)*dt_m) + psi(abs(m-1)*dt_m) - 2*psi(m*dt_m));
c0 = cov_lag(0) + 2*s_m^2;

dacf_nc = zeros(maxlag, 1);
for k = 1:maxlag
    ck = cov_lag(k);
    if k == 1, ck = ck - s_m^2; end
    dacf_nc(k) = ck / c0;
end
end

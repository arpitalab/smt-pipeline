function R = computeReturnProbability(tc, varargin)
% COMPUTERETURNPROBABILITY  Per-track "return probability" confinement metric.
%
%   R = computeReturnProbability(tc) returns a vector with one value per
%   culled track in the TrajectoryCollection tc.  Cell array of tracks is
%   also accepted in place of tc.
%
%   For each localization i in a track, define
%       N_r(i)     = #{ j ~= i : ||x_j - x_i|| <= r }
%       N_r,Ts(i)  = #{ j ~= i : ||x_j - x_i|| <= r AND |f_j - f_i|*dt <= Ts }
%   The per-localization ratio is N_r,Ts(i) / N_r(i) (set to NaN when
%   N_r(i)=0).  The per-track score R is the mean of these ratios.
%
%   Interpretation (orthogonal-ish to Pc):
%     R ~ 1   localizations within radius r are clustered in time -->
%             trajectory passes through that region once (free / drifting).
%     R << 1  spatial neighbours are spread across the whole track -->
%             particle revisits the region many times (confined / caged).
%
%   Name-value parameters:
%     'Radius'     (µm)      default 0.1     spatial neighbourhood r
%     'TimeWindow' (s)       default 0.1     temporal window Ts
%     'FrameInterval' (s)    default tc.Wrappers{1}.FrameInterval (or 0.01
%                            if a cell array was passed)
%     'MinNeighbors'         default 1       skip localizations with
%                                            N_r(i) < MinNeighbors
%     'Verbose'              default false
%
%   Example:
%     R  = computeReturnProbability(tc, 'Radius', 0.15, 'TimeWindow', 0.2);
%     Pc = tc.getPcDistribution();      % existing confinement metric
%     scatter(Pc, R); xlabel('Pc'); ylabel('Return prob');

    p = inputParser;
    addParameter(p, 'Radius',        0.1,  @(x) isscalar(x) && x > 0);
    addParameter(p, 'TimeWindow',    0.1,  @(x) isscalar(x) && x > 0);
    addParameter(p, 'FrameInterval', [],   @(x) isempty(x) || (isscalar(x) && x > 0));
    addParameter(p, 'MinNeighbors',  1,    @(x) isscalar(x) && x >= 1);
    addParameter(p, 'Verbose',       false,@islogical);
    parse(p, varargin{:});
    o = p.Results;

    if isa(tc, 'TrajectoryCollection')
        tracks = tc.getAllCulledTracks();
        if isempty(tracks)
            % Force aggregation
            tracks = {};
            for w = 1:numel(tc.Wrappers)
                tracks = [tracks, tc.Wrappers{w}.getCulledTracks()]; %#ok<AGROW>
            end
        end
        if isempty(o.FrameInterval)
            o.FrameInterval = tc.Wrappers{1}.FrameInterval;
        end
    elseif iscell(tc)
        tracks = tc;
        if isempty(o.FrameInterval)
            o.FrameInterval = 0.01;
            warning('FrameInterval not provided; assuming %.3f s.', o.FrameInterval);
        end
    else
        error('First argument must be a TrajectoryCollection or cell array of tracks.');
    end

    r2     = o.Radius^2;
    Tsfrm  = o.TimeWindow / o.FrameInterval;   % time window in frames
    nT     = numel(tracks);
    R      = nan(nT, 1);

    for k = 1:nT
        t = tracks{k};
        if size(t, 1) < 2
            continue;
        end
        xy = t(:, 1:2);
        fr = t(:, 3);
        n  = size(xy, 1);

        d2  = sum((reshape(xy, n, 1, 2) - reshape(xy, 1, n, 2)).^2, 3);
        df  = abs(fr - fr.');
        eye_n = logical(eye(n));

        within_r  = (d2 <= r2) & ~eye_n;
        within_rT = within_r & (df <= Tsfrm);

        Nr  = sum(within_r,  2);
        NrT = sum(within_rT, 2);

        valid = Nr >= o.MinNeighbors;
        if ~any(valid)
            continue;
        end
        R(k) = mean(NrT(valid) ./ Nr(valid));

        if o.Verbose && mod(k, 500) == 0
            fprintf('  track %d/%d  R=%.3f\n', k, nT, R(k));
        end
    end

    if o.Verbose
        fprintf('Return-probability: %d/%d tracks scored, median R = %.3f\n', ...
            sum(~isnan(R)), nT, median(R, 'omitnan'));
    end
end

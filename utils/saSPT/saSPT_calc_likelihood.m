function [log_L, jumps_per_track] = saSPT_calc_likelihood(tracks, D, loc_error, dt, t_exp, R)
% SASPT_CALC_LIKELIHOOD  Brownian-motion log-likelihood on (D, sigma) grid.
%
%   [log_L, jumps_per_track] = saSPT_calc_likelihood(tracks, D, loc_error, dt)
%   [...] = saSPT_calc_likelihood(tracks, D, loc_error, dt, t_exp, R)
%
%   For each subtrack, evaluates the Gaussian log-likelihood of the observed
%   step vector under the regular-Brownian-motion + isotropic-Gaussian-
%   localization-noise model (Heckert et al. 2021), at every grid point
%   (D, loc_error).  The two spatial dimensions are treated as independent,
%   so log P(dx, dy) = log P(dx) + log P(dy) and the 1D normalizer is
%   counted twice — this is the corrected version of the original
%   calc_likelihood.m, which omitted that doubling.
%
%   Motion-blur correction (Berglund 2010, Heckert 2021):
%       var(Δx)        = 2(D·dt + σ²) − 4·R·D·t_exp
%       cov(Δx, Δx_+1) = −σ² + 2·R·D·t_exp
%   with shutter coefficient R (R = 1/6 for uniform/box exposure).  When
%   t_exp is omitted (or set to 0) the static covariance is recovered.
%
%   Inputs
%     tracks    : cell array of [x, y] matrices (µm), gap-free
%     D         : 1 x LD diffusion-coefficient grid (µm² s⁻¹)
%     loc_error : 1 x LL localization-precision grid (µm)
%     dt        : sampling interval (s)
%     t_exp     : exposure time (s).  Default 0 (no blur correction).
%     R         : shutter coefficient.  Default 1/6 (uniform exposure).
%
%   Outputs
%     log_L           : LD x LL x N log-likelihood array
%     jumps_per_track : 1 x N number of one-frame jumps per subtrack

    if nargin < 5 || isempty(t_exp), t_exp = 0;   end
    if nargin < 6 || isempty(R),     R     = 1/6; end

    n_tracks = min(numel(tracks), 40000);
    LD = numel(D);
    LL = numel(loc_error);
    log_L = zeros(LD, LL, n_tracks);
    jumps_per_track = zeros(n_tracks, 1);

    LOG_2PI = log(2*pi);

    parfor itrack = 1:n_tracks
        xy = tracks{itrack};
        dx = diff(xy(:, 1));  dx(isnan(dx)) = 0;
        dy = diff(xy(:, 2));  dy(isnan(dy)) = 0;
        n  = numel(dx);
        jumps_per_track(itrack) = n;

        local = zeros(LD, LL);
        for i = 1:LD
            for j = 1:LL
                C = make_cov(D(i), loc_error(j), n, dt, t_exp, R);
                % 2D joint normalizer: −log|C| − n·log(2π)
                log_norm = 2*(n*LOG_2PI + LogDet(C));
                xll = dx' * (C \ dx);
                yll = dy' * (C \ dy);
                local(i, j) = -0.5 * (xll + yll) - 0.5 * log_norm;
            end
        end
        log_L(:, :, itrack) = local;
    end
end

function C = make_cov(D, sigma, n, dt, t_exp, R)
% Step covariance for 1D Brownian motion + Gaussian loc noise + motion
% blur under uniform exposure (Berglund 2010 / Heckert 2021).
    le2  = sigma.^2;
    blur = 2*R*D*t_exp;             % 2·R·D·t_e
    main = 2*(D*dt + le2) - 2*blur; % diagonal of step covariance
    off  = -(le2 - blur);           % lag-1 off-diagonal
    C = diag(main*ones(n,1)) ...
      + diag(off*ones(n-1,1),  1) ...
      + diag(off*ones(n-1,1), -1);
end

function y = LogDet(A)
    [U, p] = chol(A);
    if p > 0
        y = -inf;
    else
        y = 2*sum(log(diag(U)));
    end
end

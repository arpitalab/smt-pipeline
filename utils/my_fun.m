function F = my_fun(x, lags, e_msd, dt, frac)
% MY_FUN  Fractional Brownian motion MSD model for lsqnonlin fitting.
%
%   F = my_fun(x, lags, e_msd, dt, frac)
%
%   Inputs:
%     x      - parameter vector [G, sigma_sq, alpha]
%                G        : generalized diffusion coefficient (µm² / s^alpha)
%                sigma_sq : localization noise variance (µm²)
%                alpha    : anomalous exponent (0 < alpha ≤ 2)
%     lags   - integer lag indices (1, 2, 3, ...) over which to fit
%     e_msd  - empirical ensemble-mean MSD at each lag (µm²)
%     dt     - frame interval in seconds (e.g. 0.02 for 50 Hz)
%     frac   - exposure / frame_interval ratio
%                = 1     for stroboscopic illumination matching the interval
%                = e/dt  for short-pulse illumination (e.g. 0.05 for 10ms/200ms)
%
%   Output:
%     F      - vector of squared log-residuals (for lsqnonlin)
%
%   Model (Backlund et al. Phys. Rev. E):
%     b      = (|1 + frac/lag|^(2+a) + |1 - frac/lag|^(2+a) - 2) / (frac/lag)^2
%     MSD(lag) = G/((1+a)(2+a)) * ((dt*lag)^a * b - 2*(frac*dt)^a) + 2*sigma_sq
%
%   Typical usage inside getMSD():
%     fitpars = lsqnonlin(@(x) my_fun(x, (1:MaxLag)', e_msd, dt, frac), ...)

G        = x(1);
sigmasq  = x(2);
alpha    = x(3);
Delta    = 1;   % multiplier (1 for standard HILO)

b      = (abs(1 + frac ./ lags).^(2 + alpha) + ...
          abs(1 - frac ./ lags).^(2 + alpha) - 2) ./ (frac ./ lags).^2;
mlfit  = G / ((1 + alpha) * (2 + alpha)) .* ...
         ((Delta * dt .* lags).^alpha .* b - 2 * (frac * dt)^alpha) + 2 * sigmasq;

F = (log(mlfit) - log(e_msd)).^2;
end

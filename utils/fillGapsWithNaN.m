function [t_uniform, data_filled] = fillGapsWithNaN(t_nonuniform, data)
% FILLGAPSWITHNANS  Insert NaN values at missing frames in a gapped trajectory.
%
%   [t_uniform, data_filled] = fillGapsWithNaN(t_nonuniform, data)
%
%   Inputs:
%     t_nonuniform - non-uniformly spaced time vector (column or row)
%     data         - corresponding data values (same length)
%
%   Outputs:
%     t_uniform   - uniformly spaced time vector (no gaps)
%     data_filled - data with NaN inserted at missing time points

t_nonuniform = t_nonuniform(:);
data         = data(:);

if length(t_nonuniform) ~= length(data)
    error('fillGapsWithNaN: input arrays must have the same length.');
end

dt = min(diff(t_nonuniform));
if dt == 0
    % Duplicate timestamps — return as-is with a synthetic index
    t_uniform   = (1:length(t_nonuniform))';
    data_filled = data;
else
    t_uniform   = (t_nonuniform(1):dt:t_nonuniform(end))';
    data_filled = NaN(size(t_uniform));
    [~, idx_uniform, idx_nonuniform] = intersect( ...
        round(t_uniform / dt), round(t_nonuniform / dt));
    data_filled(idx_uniform) = data(idx_nonuniform);
end
end

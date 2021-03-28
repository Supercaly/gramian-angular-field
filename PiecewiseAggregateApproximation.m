function X_paa = PiecewiseAggregateApproximation(X, out_size)
% Computes Piecewise Aggregate Approximation.
%
% Parameters
% ----------
% X : array-like, shape = (n_samples, n_timestamps)
%
% output_size : Size of the returned time series. If float, it represents
%       a percentage of the size of each time series and must be
%       between 0. and 1.
%
% References
% ----------
% .. [paa](https://github.com/johannfaouzi/pyts/blob/master/pyts/approximation/paa.py)
[n_samples, n_timestamps] = size(X);

if floor(out_size) == out_size
    if out_size < 1 || out_size > n_timestamps
        throw(MException("PiecewiseAggregateApproximation:input_error", "If 'output_size' is an integer, it must be greater than or equal to 1 and lower than or equal to n_timestamps."));
    end
else
    if out_size <= 0 || out_size > 1
        throw(MException("PiecewiseAggregateApproximation:input_error", "If 'output_size' is a float, it must be greater than 0 and lower than or equal to 1 ."));
    end
    out_size = ceil(out_size * n_timestamps);
end

window_size = fix(n_timestamps / out_size);
remainder = rem(n_timestamps, out_size);
if remainder ~= 0
    window_size = window_size + 1;
end

if window_size == 1
    X_paa = X;
else
    [s,e,n_timestamps_new] = segmentation(n_timestamps, out_size);
    X_paa = zeros(n_samples, n_timestamps_new);
    for i = 1:n_samples
        for j = 1:n_timestamps_new
            % convert the index in start to matlab array space
            % the index in end is not converted because then we need
            % to do end-1 so end+1-1 makes no sense
            start_conv = s(j)+1;
            X_paa(i, j) = mean(X(i, start_conv:e(j)));
        end
    end
end
end

function [startv, endv, size] = segmentation(ts_size, n_segments)
% Compute the indices for Piecewise Agrgegate Approximation.
if floor(ts_size) ~= ts_size
    throw(MException("segmentation:input_error", "'ts_size' must be an integer."));
end
if ts_size < 2
    throw(MException("segmentation:input_error", "'ts_size' must be an integer greater than or equal to 2."));
end
if floor(n_segments) ~= n_segments
    throw(MException("segmentation:input_error", "'n_segments' must be an integer.'"));
end
if n_segments < 2
    throw(MException("segmentation:input_error", "n_segments must be greater than or equal to 2."));
end
if n_segments > ts_size
    throw(MException("segmentation:input_error", "n_segments must be lower than or equal to ts_size."));
end
 
bounds = linspace(0, ts_size, n_segments + 1);
startv = bounds(1:end-1);
endv = bounds(2:end);
size = length(startv);
end
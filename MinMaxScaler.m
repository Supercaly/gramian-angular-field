function X_new = MinMaxScaler(X, sample_range)
% Transforms samples by scaling each sample to a given range.
%
% Parameters
% ----------
% X : array-like, shape = (n_samples, n_timestamps)
%
% sample_range : tuple (min, max). Desired range of transformed data.
%
% References
% ----------
% ..[scaler](https://github.com/johannfaouzi/pyts/blob/master/pyts/preprocessing/scaler.py)
maxv = sample_range(2);
minv = sample_range(1);
X = transpose(X);
X_std = (X - min(X)) ./ (max(X) - min(X));
X_new= X_std * (maxv - minv) + minv;
X_new = transpose(X_new);
end
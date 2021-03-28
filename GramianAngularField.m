function X_new = GramianAngularField(X, image_size, method, sample_range)
% Computes the Gramian Angular Field.
% Parameters
% ----------
% X : array-like, shape = (n_samples, n_timestamps)
%
% image_size : Shape of the output images.
%       Output images are square, thus providing the size of one dimension is enough.
%
% method : 'summation' or 'difference'
%
% sample_range : Tuple (min, max)
%       Desired range of transformed data. 
%       Each sample is scaled between min and max; min must be
%       greater than or equal to -1 and max must be lower than or equal to 1.
%
% References
% ----------
% .. [gaf](https://github.com/johannfaouzi/pyts/blob/master/pyts/image/gaf.py)
[n_samples, n_timestamps] = size(X);
if image_size < 1 || image_size > n_timestamps
    throw(MException("GramianAngularField:input_error", "'image_size' must be greater than or equal to 1 and lower than or equal to 'n_timestamps'"));
end

X_paa = PiecewiseAggregateApproximation(X, image_size);
X_cos = MinMaxScaler(X_paa, sample_range);
X_sin = sqrt(clip(1 - X_cos.^2, 0, 1));
if method == "summation"
    X_new = gasf(X_cos, X_sin, n_samples, image_size);
elseif method == "difference"
    X_new = gadf(X_cos, X_sin, n_samples, image_size);
else
    throw(MException("GramianAngularField:input_error", "GramianAngularField's method must be 'summation' or 'difference'"));
end
end

function X_gasf = gasf(X_cos, X_sin, n_samples, image_size)
% Helper function that computes the step for a GASF method
X_gasf = zeros(image_size, image_size, n_samples);
for i = 1:n_samples
    cosi = X_cos(i,:);
    sini = X_sin(i,:);
    X_gasf(:,:,i)= (cosi(:)*cosi(:).') - (sini(:)*sini(:).');
end
end

function X_gadf = gadf(X_cos, X_sin, n_samples, image_size)
% Helper function that computes the step for a GADF method
X_gadf = zeros(image_size, image_size, n_samples);
for i = 1:n_samples
    cosi = X_cos(i,:);
    sini = X_sin(i,:);
    X_gadf(:,:,i)= (sini(:)*cosi(:).') - (cosi(:)*sini(:).');
end
end

function res = clip(a, a_min, a_max)
% Clip (limit) the values in an array.
%Given an interval, values outside the interval are clipped to the 
%interval edges. For example, if an interval of [0, 1] is specified, 
%values smaller than 0 become 0, and values larger than 1 become 1.
%No check is performed to ensure a_min < a_max.
%[numpy.clip](https://numpy.org/doc/stable/reference/generated/numpy.clip.html)
res = a;
res(res < a_min) = a_min;
res(res > a_max) = a_max;
end
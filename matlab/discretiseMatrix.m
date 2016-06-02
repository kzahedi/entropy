function [ m ] = discretiseMatrix(m, varargin )
%DISCRETISEMATRIX discretises a matrix
%   Required parameter:  data
%   Optional parameters: min, max, bins
%   Ouput: Vector with values in [1, bins]

p = inputParser;
p.addOptional('min',   0, @isscalar);
p.addOptional('max',   0, @isscalar);
p.addOptional('bins', 30, @isscalar);
p.parse(varargin{:});
p.Results;

min_ = p.Results.min;
max_ = p.Results.max;
bins = p.Results.bins;

if min_ == max_
    min_ = min(m);
    max_ = max(m);
end

m = min(bins, int64(1.0 + floor((m - min_) / (max_ - min_) * bins)));

end

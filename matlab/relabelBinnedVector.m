function [ r ] = relabelBinnedVector( v )
%RELABELBINNEDVECTOR Summary of this function goes here
%   Detailed explanation goes here

indices = sort(unique(v));

l = size(v, 1);
r = zeros(l, 1);

for i=1:l
    r(i,1) = find(indices == v(i,1));
end

end


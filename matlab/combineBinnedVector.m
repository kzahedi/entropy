function [ r ] = combineBinnedVector( v, bins )
%COMBINEBINNEDVECTOR Takes a vector of integers and combines them into a
%single integer
%   Detailed explanation goes here

r = 0;
for i = 1:length(v)
    r = r + v(i) * bins^(i-1);
end
end
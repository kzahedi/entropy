function [ r ] = combineBinnedMatrix( data )
%COMBINEANDRELABEL Convert a matrix into a vector
%   Call combineBinnedVector for each row of the matrix

m = max(max(data));
l = size(data,1);
r = zeros(l,1);

for i = 1:l
    r(i,1) = combineBinnedVector(data(i,:), m);
end

end
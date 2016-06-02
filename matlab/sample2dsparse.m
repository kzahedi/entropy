function [ p ] = sample2dsparse(xy)
%SAMPLE2D Emperical estimation of the joint distribution p(x,y) 
%   Detailed explanation goes here

A = double(xy);
[C,ia,ic] = unique(A,'rows');

for i=1:size(C,1)
    p(i,:)=[ C(i,:),length(find(ic==i))/size(A,1) ]; 
    % each row corresponds to one value of (x,y,z). The first two entries
    % are x,y the third entry is the number of times that it was
    % observed, divided by the total number of observations.
end

end


function [ r ] = combineAndRelabelBinnedMatrix( data )
%COMBINEANDRELABELBINNEDMATRIX Combines the rows of a discretised matrix to
%a discretised vector and relabels the entries
%   1. Make a unary random variable by combining the colums of the matrix
%   2. Beacuse the entries in the vector are used to initilaise and acces
%   matrices, we relabel the entries so that gaps in the labels

r = relabelBinnedVector(combineBinnedMatrix(data));

end


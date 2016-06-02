function [ r ] = MC_W_CMIsparse_dynamic( w2, w1, a1)
%MC_W_CMI calculates MI(W';W|A) based on the CMI function
%   Input: data, which is a matrix with 2 columns and n rows
%          column 1 = w
%          column 2 = a
%   Calculation is based on the joint distribution p(w',w,a), which is
%   handed to the CMI function

r = CMI_dynamic(w2, w1, a1);

end
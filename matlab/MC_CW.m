function [ r ] = MC_CW( w2, w1, a1 )
%MC_W_CMI calculates CIF(W -> W') - CIF(A -> W') = I(W';W) - I(W';A)
%   Input: w2, w1, a1
%   I(W';W) - I(W';A) = H(W'|A) - H(W'|W)
%   I(W';W) - I(W';A) = 
%     \sum_{w,w'} p(w',w) log_2 p(w',w) / (p(w') * p(w))
%   - \sum_{a,w'} p(,a) log_2 p(w',a) / (p(w') * p(a))
%   H(W'|A) - H(W'|W) = 
%     \sum_{a,w'} p(w',a) log_2 p(w'|a)
%   - \sum_{w,w'} p(w',w) log_2 p(w'|w)
%
% Required: p(w',a), p(a), p(w',w), p(w)

pw2w1sparse = sample2dsparse([w2 w1]);
pw2a1sparse = sample2dsparse([w2 a1]);

% pw1 = sum(pw2w1, 1);
[W1,ip,iw1] = unique(pw2w1sparse(:,2),'rows'); 
for i=1:size(W1,1)
pw1(i,:) = [W1(i,:), sum(pw2w1sparse(find(iw1==i),3))];
end

%pa1   = sum(pw2a1, 1);
[A1,ip,ia1] = unique(pw2a1sparse(:,2),'rows'); 
for i=1:size(A1,1)
pa1(i,:) = [A1(i,:), sum(pw2a1sparse(find(ia1==i),3))];
end

hw2_c_w1 = 0;
for i = 1:size(pw2w1sparse,1)
    w1 = pw2w1sparse(i,2); 
    hw2_c_w1 = hw2_c_w1 - pw2w1sparse(i,3) * log2(pw2w1sparse(i,3) / pw1(rowindex(pw1(:,1),w1),2));
end

hw2_c_a1 = 0;
for i = 1:size(pw2a1sparse,1)
    a1 = pw2a1sparse(i,2);
    hw2_c_a1 = hw2_c_a1 - pw2a1sparse(i,3) * log2(pw2a1sparse(i,3) / pa1(rowindex(pa1(:,1),a1),2));
end

r = hw2_c_a1 - hw2_c_w1;

end

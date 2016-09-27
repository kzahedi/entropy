function [ r ] = MC_CW_dynamic( w2, w1, s1, a1 )
%MC_W_CMI_dynamic calculates CIF(W -> W') - CIF(A -> W') = I(W';W) - I(W';A)
%   for each time step
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

r = zeros(length(w2), 1);

for t = 1:length(w2)
    w2t = w2(t,1);
    w1t = w1(t,1);
    a1t = a1(t,1);

    w2w1rowindex = rowindex(pw2w1sparse(:, 1:2), [w2t w1t]);
    w2a1rowindex = rowindex(pw2a1sparse(:, 1:2), [w2t a1t]);
    w1rowindex   = rowindex(pw1(:, 1), w1t);  
    a1rowindex   = rowindex(pa1(:, 1), a1t);

    hw2_c_w1     = -log2(pw2w1sparse(w2w1rowindex, 3) / (pw1(w1rowindex, 2)));
    hw2_c_a1     = -log2(pw2a1sparse(w2a1rowindex, 3) / (pa1(a1rowindex, 2)));

    r(t,1) = hw2_c_a1 - hw2_c_w1;
end

end

function [ r ] = MC_MI( w2, w1, s1, a1 )
%MC_MI Calculates I(W';W) - I(A;S) as sum of entropies
%   I(W';W) - I(A;S) = H(W') - H(W'|A) + H(A) - H(A|S)
%   H(W')       = \sum p(w') log2( p(w') )
%   H(W'|W)     = \sum p(w',w) log2( p(w'|w) )
%   H(A)        = \sum p(a) log2( p(a) )
%   H(A|S)      = \sum p(a,s) log2( p(a|s) )
%   p(a,s)      = by sampling
%   p(w',w)     = by sampling
%   p(w')       = \sum_w p(w')
%   p(a)        = \sum_s p(a,s)

% pw2w1a1s1 = sample4d([w2 w1 a1 s1]);

pw2w1sparse = sample2dsparse([w2 w1]);
pa1s1sparse = sample2dsparse([a1 s1]);

% pw1 = sum(pw2w1, 1);
[W1,ip,iw1] = unique(pw2w1sparse(:,2),'rows'); 
for i=1:size(W1,1)
pw1(i,:) = [W1(i,:), sum(pw2w1sparse(find(iw1==i),3))];
end

% pw2   = sum(pw2w1, 2);
[W2,ip,iw2] = unique(pw2w1sparse(:,1),'rows'); 
for i=1:size(W2,1)
pw2(i,:) = [W2(i,:), sum(pw2w1sparse(find(iw2==i),3))];
end

% pa1   = sum(pa1s1, 2);
[A1,ip,ia1] = unique(pa1s1sparse(:,1),'rows'); 
for i=1:size(A1,1)
pa1(i,:) = [A1(i,:), sum(pa1s1sparse(find(ia1==i),3))];
end

%ps1   = sum(pa1s1, 1);
[S1,ip,is1] = unique(pa1s1sparse(:,2),'rows'); 
for i=1:size(S1,1)
ps1(i,:) = [S1(i,:), sum(pa1s1sparse(find(is1==i),3))];
end

hw2 = 0;
for i = 1:size(pw2,1)
        hw2 = hw2 - pw2(i,2) * log2(pw2(i,2));
end

ha1 = 0;
for i = 1:size(pa1,1)
        ha1 = ha1 - pa1(i,2) * log2(pa1(i,2));
end

hw2_c_w1 = 0;
for i = 1:size(pw2w1sparse,1)
    w2 = pw2w1sparse(i,1);
    w1 = pw2w1sparse(i,2); 
        hw2_c_w1 = hw2_c_w1 - pw2w1sparse(i,3) * log2(pw2w1sparse(i,3) / pw1(rowindex(pw1(:,1),w1),2));
end

ha1_c_s1 = 0;
for i = 1:size(pa1s1sparse,1)
    a1 = pa1s1sparse(i,1);
    s1 = pa1s1sparse(i,2);
            ha1_c_s1 = ha1_c_s1 - pa1s1sparse(i,3) * log2(pa1s1sparse(i,3) / ps1(rowindex(ps1(:,1),s1),2));
end

r = hw2 - hw2_c_w1 - ha1 + ha1_c_s1;
end
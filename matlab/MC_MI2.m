function [ r ] = MC_MI2( w2, w1, s1, a1 )
%MC_MI1 I(W';W) - I(A;S)
%   I(W';W) - I(A;S) = H(W') - H(W'|A) + H(A) - H(A|S)
%   H(W')       = \sum p(w') log2( p(w') )
%   H(W'|W)     = \sum p(w',w) log2( p(w'|w) )
%   H(A)        = \sum p(a) log2( p(a) )
%   H(A|S)      = \sum p(a,s) log2( p(a|s) )
%   p(a,s)      = by sampling
%   p(w',w)     = by sampling
%   p(w')       = \sum_w p(w')
%   p(a)        = \sum_s p(a,s)

w2max = max(w2);
w1max = max(w1);
a1max = max(a1);
s1max = max(s1);

pw2w1    = zeros(w2max, w1max);
pw2      = zeros(w2max, 1);
pw1      = zeros(w1max, 1);
pw2_c_w1 = zeros(w2max, w1max);

pa1s1    = zeros(a1max, s1max);
pa1      = zeros(a1max, 1);
ps1      = zeros(s1max, 1);
pa1_c_s1 = zeros(a1max, s1max);

for t = 1:length(w2)
    w2index = w2(t);
    w1index = w1(t);
    a1index = a1(t);
    s1index = s1(t);

    pw2w1(w2index, w1index) = pw2w1(w2index, w1index) + 1.0;
    pw2(w2index)            = pw2(w2index)            + 1.0;
    pw1(w1index)            = pw1(w1index)            + 1.0;
    pa1s1(a1index, s1index) = pa1s1(a1index, s1index) + 1.0;
    pa1(a1index)            = pa1(a1index)            + 1.0;
    ps1(s1index)            = ps1(s1index)            + 1.0;

end

pw2w1 = pw2w1 / length(w2);
pw2   = pw2   / length(w2);
pw1   = pw1   / length(w2);

pa1s1 = pa1s1 / length(w2);
pa1   = pa1   / length(w2);
ps1   = ps1   / length(w2);

for w2index = 1:w2max
    for w1index = 1:w1max
        pw2_c_w1(w2index, w1index) = pw2w1(w2index, w1index) / pw1(w1index);
    end
end

for a1index = 1:a1max
    for s1index = 1:s1max
        pa1_c_s1(a1index, s1index) = pa1s1(a1index, s1index) / ps1(s1index);
    end
end

h_w2 = 0;
for w2index = 1:w2max
    if pw2(w2index) > 0.0
        h_w2 = h_w2 - pw2(w2index) * log2( pw2(w2index) );
    end
end

h_a1 = 0;
for a1index = 1:a1max
    if pa1(a1index) > 0.0
        h_a1 = h_a1 - pa1(a1index) * log2( pa1(a1index) );
    end
end

h_w2_c_w1 = 0;
for w2index = 1:w2max
    for w1index = 1:w1max
        if pw2w1(w2index, w1index) > 0.0 && pw2_c_w1(w2index, w1index) > 0.0
            h_w2_c_w1 = h_w2_c_w1 - pw2w1(w2index, w1index) * log2( pw2_c_w1(w2index, w1index));
        end
    end
end

h_a1_c_s1 = 0;
for a1index = 1:a1max
    for s1index = 1:s1max
        if pa1s1(a1index, s1index) > 0.0 && pa1_c_s1(a1index, s1index) > 0.0
            h_a1_c_s1 = h_a1_c_s1 - pa1s1(a1index, s1index) * log2( pa1_c_s1(a1index, s1index));
        end
    end
end


i_w2_w1 = h_w2 - h_w2_c_w1;

i_a1_s1 = h_a1 - h_a1_c_s1;

r = i_w2_w1 - i_a1_s1;

end


function [ r ] = MC_MI1( w2, w1, s1, a1 )
%MC_MI1 I(W';W) - I(A;S)
%   I(W';W) = \sum_{w', w} p(w',w) log2 (p(w',w) / (p(w') p(w))
%   I(A;S)  = \sum_{a, s}  p(a,s)  log2 (p(a,s)  / (p(a)  p(s))

w2max = max(w2);
w1max = max(w1);
a1max = max(a1);
s1max = max(s1);

pw2w1 = zeros(w2max, w1max);
pw2   = zeros(w2max, 1);
pw1   = zeros(w1max, 1);

pa1s1 = zeros(a1max, s1max);
pa1   = zeros(a1max, 1);
ps1   = zeros(s1max, 1);

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


i_w2_w1 = 0;
for w2index = 1:w2max
    for w1index = 1:w1max
        if pw2w1(w2index, w1index) > 0.0 && ...
           pw2(w2index)            > 0.0 && ...
           pw1(w1index)            > 0.0
       i_w2_w1 = i_w2_w1 + pw2w1(w2index, w1index) * log2( pw2w1(w2index, w1index) / (pw2(w2index) * pw1(w1index)));
        end
    end
end

i_a1_s1 = 0;
for a1index = 1:a1max
    for s1index = 1:s1max
        if pa1s1(a1index, s1index) > 0.0 && ...
           pa1(a1index)            > 0.0 && ...
           ps1(s1index)            > 0.0
       i_a1_s1 = i_a1_s1 + pa1s1(a1index, s1index) * log2( pa1s1(a1index, s1index) / (pa1(a1index) * ps1(s1index)));
        end
    end
end

r = i_w2_w1 - i_a1_s1;

end


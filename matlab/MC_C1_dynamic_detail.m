function [ r ] = MC_C1_dynamic_detail( w2, w1, s1, a1 )
%MC_C = CIF(W -> W') - CIF(A -> W')
% CIF(W -> W') = p(w',w) log2( p(w'|w) / p(w') ) = p(w',w) log2( p(w',w) / (p(w') p(w)) )
% CIF(A -> W') = \sum_{w',a} p(w'|do(a)) p(a) log2(p(w'|do(a)) / \hat{p}(w))
% \hat{p}(w')  = \sum_a p(w'|do(a))p(a) % post-interventionale verteilung von a
% p(w'|do(a))  = \sum_s p(w'|s,a) p(s)
% p(w'|s,a)    = p(w',s,a) / p(s,a)

w2max = max(w2);
w1max = max(w1);
a1max = max(a1);
s1max = max(s1);

pw2s1a1 = zeros(w2max, s1max, a1max);
ps1a1   = zeros(s1max, a1max);
ps1     = zeros(s1max, 1);
pa1     = zeros(a1max, 1);
pw2w1   = zeros(w2max, w1max);
pw2     = zeros(w2max, 1);
pw1     = zeros(w1max, 1);
hatpw2  = zeros(w2max, 1);

cifww_values  = zeros(w2max, w1max, a1max);
cifwa_values  = zeros(w2max, w1max, a1max);

pw2_c_s1a1  = zeros(w2max, s1max, a1max);
pw2_c_do_a1 = zeros(w2max, a1max);

for t = 1:length(w2)
    w2index = w2(t);
    w1index = w1(t);
    a1index = a1(t);
    s1index = s1(t);
    
    pw2s1a1(w2index, s1index, a1index) = pw2s1a1(w2index, s1index, a1index) + 1.0;
    ps1a1(s1index, a1index)            = ps1a1(s1index, a1index)            + 1.0;
    ps1(s1index,1)                     = ps1(s1index,1)                     + 1.0;
    pa1(a1index,1)                     = pa1(a1index,1)                     + 1.0;
    pw2w1(w2index, w1index)            = pw2w1(w2index, w1index)            + 1.0;
    pw2(w2index,1)                     = pw2(w2index,1)                     + 1.0;
    pw1(w1index,1)                     = pw1(w1index,1)                     + 1.0;
end

pw2s1a1 = pw2s1a1 / length(w2);
ps1a1   = ps1a1   / length(w2);
ps1     = ps1     / length(w2);
pa1     = pa1     / length(w2);
pw2w1   = pw2w1   / length(w2);
pw2     = pw2     / length(w2);
pw1     = pw1     / length(w2);

for w2index = 1:w2max
    for a1index = 1:a1max
        for s1index = 1:s1max
            if ps1a1(s1index, a1index) > 0.0
                pw2_c_s1a1(w2index, s1index, a1index) = pw2s1a1(w2index, s1index, a1index) / ps1a1(s1index, a1index);
                if isnan(pw2_c_s1a1(w2index, s1index, a1index))
                    fprintf('pw2_c_s1a1(%d, %d, %d) is nan\n', w2index, s1index, a1index);
                end
                if isnan(pw2s1a1(w2index, s1index, a1index))
                    fprintf('pw2s1a1(%d, %d, %d) is nan\n', w2index, s1index, a1index);
                end
                if isnan(ps1a1(s1index, a1index))
                    fprintf('ps1a1(%d, %d) is nan\n', s1index, a1index);
                end
            end
        end
    end
end

% p(w'|do(a))  = \sum_s p(w'|s,a) p(s)
for w2index = 1:w2max
    for a1index = 1:a1max
        for s1index = 1:s1max
            pw2_c_do_a1(w2index, a1index) = pw2_c_do_a1(w2index, a1index) + pw2_c_s1a1(w2index, s1index, a1index) * ps1(s1index,1);
        end
    end
end

% \hat{p}(w')  = \sum_a p(w'|do(a))p(a) % post-interventionale verteilung von a
for w2index = 1:w2max
    for a1index = 1:a1max
        hatpw2(w2index,1) = hatpw2(w2index,1) + pw2_c_do_a1(w2index, a1index) * pa1(a1index,1);
    end
end

% CIF(W -> W') = p(w',w) log2( p(w'|w) / p(w') ) = p(w',w) log2( p(w',w) / (p(w') p(w)) )
% CIF(A -> W') = \sum_{w',a} p(w'|do(a)) p(a) log2(p(w'|do(a)) / \hat{p}(w'))


% cif_w2_to_w1 = 0;
% for w2index = 1:w2max
%     for w1index = 1:w1max
%         if pw2(w2index) > 0.0 && pw1(w1index) > 0.0 && pw2w1(w2index, w1index) > 0.0
%             cif_w2_to_w1 = cif_w2_to_w1 + pw2w1(w2index, w1index) * log2( pw2w1(w2index, w1index) / (pw2(w2index,1) * pw1(w1index,1)));
%         end
%     end
% end

for w2index = 1:w2max
    for w1index = 1:w1max
        for a1index = 1:a1max
            if hatpw2(w2index)                    > 0.0 && ...
                    pw2_c_do_a1(w2index, a1index) > 0.0 && ...
                    pw2(w2index)                  > 0.0 && ...
                    pw1(w1index)                  > 0.0 && ...
                    pw2w1(w2index, w1index)       > 0.0
                
                cif_w2_to_w1 = log2( pw2w1(w2index, w1index) / (pw2(w2index,1) * pw1(w1index,1)));
                cif_w2_to_a1 = log2( pw2_c_do_a1(w2index, a1index) / (hatpw2(w2index,1)));
                
                cifww_values(w2index, w1index, a1index) = cif_w2_to_w1;
                cifwa_values(w2index, w1index, a1index) = cif_w2_to_a1;
            end
        end
    end
end

r = zeros(length(w2), 3);

for t = 1:length(w2)
    w2t = w2(t);
    w1t = w1(t);
    a1t = a1(t);
    r(t, 1) = cifww_values(w2index, w1index, a1index);
    r(t, 2) = cifwa_values(w2index, w1index, a1index);
    r(t, 3) = cifww_values(w2index, w1index, a1index) - cifwa_values(w2index, w1index, a1index);
end

end
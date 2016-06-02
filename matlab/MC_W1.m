function [ r ] = MC_W1( w2, w1, a1 )
%MC_W_CMI calculates MI(W';W|A)
%   Input: w2, w1, a1
%   I(W';W|A) = \sum_{w', w, a) p(w', w, a) log_2 (p(w'|w,a) / p(w'|a))
%   p(w'|w,a) = p(w', w, a) / p(w, a)
%   p(w'|a)   = p(w', a)    / p(a)

w2max = max(w2);
w1max = max(w1);
a1max = max(a1);

pw2w1a1    = zeros(w2max, w1max, a1max); % p(w', w, a)
pw1a1      = zeros(w1max, a1max);        % p(w, a)
pw2_c_w1a1 = zeros(w2max, w1max, a1max); % p(w'| w, a)
pw2a1      = zeros(w2max, a1max);        % p(w', a)
pa1        = zeros(a1max, 1);            % p(a)
pw2_c_a1   = zeros(w2max, a1max);        % p(w'|a)

%
% p(w', w, a)
% p(w, a)
% p(w', a)
% p(a)
%
for t = 1:length(w1)
    w2index = w2(t);
    w1index = w1(t);
    a1index = a1(t);
    pw2w1a1(w2index, w1index, a1index) = pw2w1a1(w2index, w1index, a1index) + 1.0;
    pw1a1(w1index, a1index)            = pw1a1(w1index, a1index) + 1.0;
    pw2a1(w2index, a1index)            = pw2a1(w2index, a1index) + 1.0;
    pa1(a1index,1)                     = pa1(a1index) + 1.0;
end
pw2w1a1 = pw2w1a1 / length(w1);
pw1a1   = pw1a1   / length(w1);
pw2a1   = pw2a1   / length(w1);
pa1     = pa1     / length(w1);

%
% p(w'|w,a) = p(w',w,a) / p(w, a)
%
for w2index = 1:w2max
    for w1index = 1:w1max
        for a1index = 1:a1max
          if pw1a1(w1index, a1index) > 0.0
            pw2_c_w1a1(w2index, w1index, a1index) = pw2w1a1(w2index, w1index, a1index) / pw1a1(w1index, a1index);
          end
        end
    end
end

%
% p(w'|a) = p(w',a) / p(a)
%
for w2index = 1:w2max
    for a1index = 1:a1max
      if pa1(a1index) > 0.0
        pw2_c_a1(w2index, a1index) = pw2a1(w2index, a1index) / pa1(a1index,1);
      end
    end
    
end

r = 0;
for w2index = 1:w2max
    for w1index = 1:w1max
        for a1index = 1:a1max
            if pw2w1a1(w2index, w1index, a1index)    > 0.0 && ...
               pw2_c_w1a1(w2index, w1index, a1index) > 0.0 && ...
               pw2_c_a1(w2index, a1index)            > 0.0
           r = r + pw2w1a1(w2index, w1index, a1index) * round(log2( pw2_c_w1a1(w2index, w1index, a1index) / pw2_c_a1(w2index, a1index)),4);
            end
        end
    end
end

end

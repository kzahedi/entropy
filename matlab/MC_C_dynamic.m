function [ r ] = MC_C_dynamic( w2, w1, s1, a1 )
%MC_C = CIF(W -> W') - CIF(A -> W')
% CIF(W -> W') = p(w',w) log2( p(w'|w) / p(w') ) = p(w',w) log2( p(w',w) / (p(w') p(w)) )
% CIF(A -> W') = \sum_{w',a} p(w'|do(a)) p(a) log2(p(w'|do(a)) / \hat{p}(w))
% \hat{p}(W)   = \sum_a p(w'|do(a))p(a) % post-interventionale verteilung von a
% p(w'|do(a))  = \sum_s p(w'|s,a) p(s)
% p(w'|s,a)    = p(w',s,a) / p(s,a)


pw2w1 = sample2dsparse([w2 w1]);
pw2a1 = sample2dsparse([w2 a1]);

% pw1  = sum(pw2w1, 1);
[W1,ip,iw1] = unique(pw2w1(:,2),'rows');
for i=1:size(W1,1)
    pw1(i,:) = [W1(i,:), sum(pw2w1(find(iw1==i),3))];
end

% pw2  = sum(pw2w1, 2);
[W2,ip,iw2] = unique(pw2w1(:,1),'rows');
for i=1:size(W2,1)
    pw2(i,:) = [W2(i,:), sum(pw2w1(find(iw2==i),3))];
end

% pa1  = sum(pw2a1, 2);
[A1,ip,ia1] = unique(pw2a1(:,2),'rows');
for i=1:size(A1,1)
    pa1(i,:) = [A1(i,:), sum(pw2a1(find(ia1==i),3))];
end


pw2s1a1 = sample3dsparse([w2 s1 a1]);

% ps1  = sum(sum(pw2s1a1,1),3);
[S1,ip,is1] = unique(pw2s1a1(:,2),'rows');
for i=1:size(S1,1)
    ps1(i,:) = [S1(i,:), sum(pw2s1a1(find(is1==i),4))];
end

ps1a1 = sample2dsparse([s1 a1]);

% p(w'|do(a))  = \sum_s p(w'|s,a) p(s) = \sum_s (p(w',s,a) / p(s,a)) p(s)
for i=1:size(pw2a1,1)
    w2index = pw2a1(i,1);
    a1index = pw2a1(i,2);
    pw2S1a1 = pw2s1a1(rowindex(pw2s1a1(:,[1 3]),[w2index a1index]), [2 4]);
    pS1a1   = ps1a1(rowindex(ps1a1(:,2),a1index),[1 3]);
    
    val=0;
    for j=1:size(ps1,1)
        s1index = ps1(j,1);
        if ~isempty(rowindex(pw2S1a1(:,1),s1index))
            val = val + (pw2S1a1(rowindex(pw2S1a1(:,1),s1index),2) / pS1a1(rowindex(pS1a1(:,1),s1index),2) ) * ps1(j,2);
        end
    end
    pw2doa(i,:) = [w2index, a1index, val];
end


% calculate the measure for each time step t

v = [];

for t = 1:length(w1)
    w2t = w2(t);
    w1t = w1(t);
    a1t = a1(t);
    
    w1index   = rowindex(pw1(:,1), w1t);
    w2index   = rowindex(pw2(:,1), w2t);
    w2w1index = rowindex(pw2w1(:,  [1 2]), [w2t, w1t]);
    w2a1index = rowindex(pw2doa(:, [1 2]), [w2t, a1t]);
    
    hatpw = 0;
    for j=1:size(pa1,1)
        a1tmp=pa1(j,1);
        if ~isempty(rowindex(pw2doa(:,[1 2]),[w2t a1tmp]))
            hatpw = hatpw + pw2doa(rowindex(pw2doa(:,[1 2]),[w2t a1tmp]),3) * pa1(j,2);
        end
    end
    
    cifw2w1 = log2( pw2w1(w2w1index,3)  / (pw2(w2index,2) * pw1(w1index,2)));
    cifw2a1 = log2( pw2doa(w2a1index,3) / hatpw );
    
    nv = cifw2w1 - cifw2a1;
    found = 0;
    for x = 1:size(v,1)
        if sum(v(x,1:3) == [w2t, w1t, a1t]) == 3
            found = 1;
            nv = v(x, 4);
        end
    end

    if found == 0
        v(t, :)= [double(w2t), double(w1t), double(a1t), cifw2w1 - cifw2a1];
    end
    
    r(t,1) = nv;

end

end
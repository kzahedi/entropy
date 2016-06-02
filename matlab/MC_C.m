function [ r ] = MC_C( w2, w1, s1, a1 )
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
    w2 = pw2a1(i,1);
    a1 = pw2a1(i,2);
    pw2S1a1 = pw2s1a1(rowindex(pw2s1a1(:,[1 3]),[w2 a1]), [2 4]); 
    pS1a1   = ps1a1(rowindex(ps1a1(:,2),a1),[1 3]);

val=0;    
    for j=1:size(ps1,1)
        s1 = ps1(j,1);
        if ~isempty(rowindex(pw2S1a1(:,1),s1))            
val = val + (pw2S1a1(rowindex(pw2S1a1(:,1),s1),2) / pS1a1(rowindex(pS1a1(:,1),s1),2) ) * ps1(j,2);
        end
    end
    pw2doa(i,:) = [w2,a1, val];
end


% CIF(A -> W') = \sum_{w',a} p(w'|do(a)) p(a) log2(p(w'|do(a)) / \hat{p}(w))
cifw2a1 = 0;

for i = 1:size(pw2doa,1)
w2 = pw2doa(i,1);

    hatpw = 0;
    for j=1:size(pa1,1)
        a1=pa1(j,1);
        if ~isempty(rowindex(pw2doa(:,[1 2]),[w2 a1]))
        hatpw = hatpw + pw2doa(rowindex(pw2doa(:,[1 2]),[w2 a1]),3) * pa1(j,2);
        end
    end
    
    a1 = pw2doa(i,2);
    
    cifw2a1 = cifw2a1 ...
                    + pw2doa(i,3) * pa1(rowindex(pa1(:,1),a1),2) ...
                    * ( log2( pw2doa(i,3) / hatpw ));
end


% CIF(W -> W') = p(w',w) log2( p(w'|w) / p(w') ) = p(w',w) log2( p(w',w) / (p(w') p(w)) )
cifw2w1 = 0;

for i = 1:size(pw2w1,1)
    w2 = pw2w1(i,1);
    w1 = pw2w1(i,2);
    cifw2w1 = cifw2w1 ...
                    + pw2w1(i,3) ...
                    * ( log2( pw2w1(i,3) / (pw2(rowindex(pw2(:,1),w2),2) * pw1(rowindex(pw1(:,1),w1),2)) ));
end




r = cifw2w1 - cifw2a1;

end
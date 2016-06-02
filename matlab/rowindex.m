function I = rowindex(M,m)
        % return the set of indices i with M(i,:)=m
        I=1:size(M,1);
        for s=1:length(m)
            I = intersect(I,find(M(:,s)==m(s)));
        end
end

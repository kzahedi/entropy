function [ r ] = CMIsparse( pxyzsparse )
%CMIsparse Conditional Mutual Information MI(X;Y|Z)
%   Input is a matrix M(x,y,z,p(x,y,z))
%   Calculation is:
%   MI(X;Y|Z) = \sum p(x,y,z) log (p(x,y|z) / (p(x|z) * p(y|z)))
%   p(x,y|z)  = p(x,y,z) / p(z)
%   p(x|z)    = sum_y p(x,y,z) / p(z)
%   p(y|z)    = sum_x p(x,y,z) / p(z)
%   p(z)      = sum_{x,y} p(x,y,z)

r = 0;

% pz  = sum(sum(pxyz, 2), 1);
[Z,ip,iz] = unique(pxyzsparse(:,3),'rows'); 
for i=1:size(Z,1)
pz(i,:) = [Z(i,:), sum(pxyzsparse(find(iz==i),4))];
end

% pxz = sum(pxyz, 2);
[XZ,ip,ixz] = unique(pxyzsparse(:,[1 3]),'rows'); 
for i=1:size(XZ,1)
pxz(i,:) = [XZ(i,:), sum(pxyzsparse(find(ixz==i),4))];
end

% pyz = sum(pxyz, 1);
[YZ,ip,iyz] = unique(pxyzsparse(:,[2 3]),'rows'); 
for i=1:size(YZ,1)
pyz(i,:) = [YZ(i,:), sum(pxyzsparse(find(iyz==i),4))];
end


for i = 1:size(pxyzsparse,1)
    x = pxyzsparse(i,1);
    y = pxyzsparse(i,2);
    z = pxyzsparse(i,3);
    r = r ...
      + pxyzsparse(i,4) ...
      * ( log2( pxyzsparse(i,4) / pz(rowindex(pz(:,1),z),2) )...
      - log2( (pxz(rowindex(pxz(:,1:2),[x z]),3) / pz(rowindex(pz(:,1),z),2) ) ...
      * ( pyz(rowindex(pyz(:,1:2),[y z]),3) / pz(rowindex(pz(:,1),z),2) ) ));
end


end

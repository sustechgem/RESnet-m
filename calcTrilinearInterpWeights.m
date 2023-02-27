function weights = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,points) 
% A function in the package "RESnet-m" 
% Core function of the trilinear interpolation. Calculate the weights of 
% the eight neighboring nodes surrounding a given point in a lattice grid.
%
% function weights = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,points)
% INPUT
%     nodeX,nodeY,nodeZ: node locations in X, Y, Z of a rectilinear mesh
%     points: a Npoints x 3 matrix specifying the X, Y, Z coordinates of points
% OUTPUT
%     weights: a Nnodes x Npoints sparse matrix of the calculated weights
% NOTE
%     weights: distribute the point's contrbution to neighboring grid nodes
%     weights' (transpose): estimate the point's value from the neighboring
%         nodes

% Prep
Nnx = length(nodeX);
Nny = length(nodeY);
Nnz = length(nodeZ);
Npoint = size(points,1);
x = points(:,1);
y = points(:,2);
z = points(:,3);

% Snap out-of-region points to the nearest boundary
x = max(x,nodeX(1));
x = min(x,nodeX(end));
y = max(y,nodeY(1));
y = min(y,nodeY(end));
z = min(z,nodeZ(1));
z = max(z,nodeZ(end));

% Nodes-to-points distance in x, y, z-direction
nxd = abs(reshape(nodeX,1,[]) - x);
nyd = abs(reshape(nodeY,1,[]) - y);
nzd = abs(reshape(nodeZ,1,[]) - z);

% Find the two nearest nodes and the distances to the nodes in x, y, z-direction
[xd, xn] = mink(nxd,2,2);
ind = xn(:,1)>xn(:,2);
tmp1 = xn(ind,1); tmp2 = xn(ind,2); xn(ind,1) = tmp2; xn(ind,2) = tmp1;
tmp1 = xd(ind,1); tmp2 = xd(ind,2); xd(ind,1) = tmp2; xd(ind,2) = tmp1;

[yd, yn] = mink(nyd,2,2);
ind = yn(:,1)>yn(:,2);
tmp1 = yn(ind,1); tmp2 = yn(ind,2); yn(ind,1) = tmp2; yn(ind,2) = tmp1;
tmp1 = yd(ind,1); tmp2 = yd(ind,2); yd(ind,1) = tmp2; yd(ind,2) = tmp1;

[zd, zn] = mink(nzd,2,2);
ind = zn(:,1)>zn(:,2);
tmp1 = zn(ind,1); tmp2 = zn(ind,2); zn(ind,1) = tmp2; zn(ind,2) = tmp1;
tmp1 = zd(ind,1); tmp2 = zd(ind,2); zd(ind,1) = tmp2; zd(ind,2) = tmp1;

% Calculate weights of eight neighboring nodes
sxd = sum(xd,2);
syd = sum(yd,2);
szd = sum(zd,2);
w1 = xd(:,2)./sxd .* yd(:,2)./syd .* zd(:,2)./szd; 
w2 = xd(:,2)./sxd .* yd(:,2)./syd .* zd(:,1)./szd;
w3 = xd(:,1)./sxd .* yd(:,2)./syd .* zd(:,2)./szd;
w4 = xd(:,1)./sxd .* yd(:,2)./syd .* zd(:,1)./szd;
w5 = xd(:,2)./sxd .* yd(:,1)./syd .* zd(:,2)./szd; 
w6 = xd(:,2)./sxd .* yd(:,1)./syd .* zd(:,1)./szd;
w7 = xd(:,1)./sxd .* yd(:,1)./syd .* zd(:,2)./szd;
w8 = xd(:,1)./sxd .* yd(:,1)./syd .* zd(:,1)./szd;

% Calculate indices of eight neighboring nodes
n1 = (Nnx*Nnz) * (yn(:,1)-1) + Nnz * (xn(:,1)-1) + zn(:,1);
n2 = (Nnx*Nnz) * (yn(:,1)-1) + Nnz * (xn(:,1)-1) + zn(:,2);
n3 = (Nnx*Nnz) * (yn(:,1)-1) + Nnz * (xn(:,2)-1) + zn(:,1);
n4 = (Nnx*Nnz) * (yn(:,1)-1) + Nnz * (xn(:,2)-1) + zn(:,2);
n5 = (Nnx*Nnz) * (yn(:,2)-1) + Nnz * (xn(:,1)-1) + zn(:,1);
n6 = (Nnx*Nnz) * (yn(:,2)-1) + Nnz * (xn(:,1)-1) + zn(:,2);
n7 = (Nnx*Nnz) * (yn(:,2)-1) + Nnz * (xn(:,2)-1) + zn(:,1);
n8 = (Nnx*Nnz) * (yn(:,2)-1) + Nnz * (xn(:,2)-1) + zn(:,2);

% Put weights in a sparse matrix: Nnodes x Npoints
s1 = sparse(n1,1:Npoint,w1,Nnx*Nny*Nnz,Npoint);
s2 = sparse(n2,1:Npoint,w2,Nnx*Nny*Nnz,Npoint);
s3 = sparse(n3,1:Npoint,w3,Nnx*Nny*Nnz,Npoint);
s4 = sparse(n4,1:Npoint,w4,Nnx*Nny*Nnz,Npoint);
s5 = sparse(n5,1:Npoint,w5,Nnx*Nny*Nnz,Npoint);
s6 = sparse(n6,1:Npoint,w6,Nnx*Nny*Nnz,Npoint);
s7 = sparse(n7,1:Npoint,w7,Nnx*Nny*Nnz,Npoint);
s8 = sparse(n8,1:Npoint,w8,Nnx*Nny*Nnz,Npoint);

% Assemble
weights = s1 + s2 + s3 + s4 + s5 + s6 + s7 + s8;


end
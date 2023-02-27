function [cellVal,faceVal,edgeVal] = ...
    makeRectMeshModelBlocks(nodeX,nodeY,nodeZ,blkLoc,blkVal,bkgCellVal,bkgFaceVal,bkgEdgeVal)
% A function in the package "RESnet-m" 
% Make edgeCon, faceCon, cellCon models using blocks
%
% function [cellVal,faceVal,edgeVal] = ...
%     makeRectMeshModelBlocks(nodeX,nodeY,nodeZ,blkLoc,blkVal,bkgCellVal,bkgFaceVal,bkgEdgeVal)
% INPUT
%     nodeX,nodeY,nodeZ: node locations in X, Y, Z of a rectilinear mesh
%     blkLoc: a Nblock x 6 matrix whose columns are [xmin xmax ymin ymax
%         zmax zmin] specifying the range of the block; if the range of any 
%         dimension is zero, that dimension vanishes to represent a thin 
%         object; one dimension vanishes for 2D sheet object; two dimensions
%         vanish for 1D line object; point object not allowed.
%     blkVal: a vector specifying the physical property values of the blocks.
%     bkgCellVal: background model values defined at cell centers; for the 
%         description of existing model structure; can be a scalar for the 
%         whole-space or a vector; treat empty [] as 0.
%     bkgFaceVal: background model values defined at cell faces; for the 
%         description of existing model structure; can be a scalar for the 
%         whole-space or a vector; treat empty [] as 0.
%     bkgEdgeVal: background model values defined at cell edges; for the 
%         description of existing model structure; can be a scalar for the 
%         whole-space or a vector; treat empty [] as 0.
% OUTPUT
%     cellVal: a vector of physical property values defined on all cells
%     faceVal: a vector of physical property values defined on all faces (cellVal * thickness)
%     edgeVal: a vector of physical property values defined on all edges (cellVal * cross-sectional area)
% NOTE
%     Entries in edgeCon and faceCon are directional. 
%     First level ordering: for directional objects (edge and face's normal), 
%         follow x, y, z-orientation ordering
%     Second level ordering: for non-directional objects (cell) and within
%         a particular orientation, count in z (top to bottom), then x 
%         (left to right), then y (front to back)

% Prep
Nx = length(nodeX);
Ny = length(nodeY);
Nz = length(nodeZ);
Nedges = (Nx-1) * Ny * Nz + Nx * (Ny-1) * Nz + Nx * Ny * (Nz-1);
Nfaces = (Nx-1) * Ny * (Nz-1) + Nx * (Ny-1) * (Nz-1) + (Nx-1) * (Ny-1) * Nz;
Ncells = (Nx-1) * (Ny-1) * (Nz-1);

if isempty(bkgCellVal)
    bkgCellVal = 0;
end
if isempty(bkgFaceVal)
    bkgFaceVal = 0;
end
if isempty(bkgEdgeVal)
    bkgEdgeVal = 0;
end
edgeVal = zeros(Nedges,1) + bkgEdgeVal;
faceVal = zeros(Nfaces,1) + bkgFaceVal;
cellVal = zeros(Ncells,1) + bkgCellVal;

% Get connectivity lists from mesh definition
[nodes, edges, ~, faces, ~, cells, ~] = ...
                                formRectMeshConnectivity(nodeX,nodeY,nodeZ);

% Replace inf with the outmost boundary
blkLoc(blkLoc(:,1)==-inf,1) = nodeX(1);
blkLoc(blkLoc(:,2)==inf,2) = nodeX(end);
blkLoc(blkLoc(:,3)==-inf,3) = nodeY(1);
blkLoc(blkLoc(:,4)==inf,4) = nodeY(end);
blkLoc(blkLoc(:,5)==inf,5) = nodeZ(1);
blkLoc(blkLoc(:,6)==-inf,6) = nodeZ(end);

% Pre-screen to identify object types: 1D, 2D or 3D
dim = ~[abs(blkLoc(:,1)-blkLoc(:,2))==0 ...
      abs(blkLoc(:,3)-blkLoc(:,4))==0 ...
      abs(blkLoc(:,5)-blkLoc(:,6))==0  ]; % 0 indicates that dimension vanished
objType = sum(dim,2); % dimensionality = 3 for volume, 2 for sheet, 1 for string

% Object center positions
edgesCenter = 1/2 * (nodes(edges(:,1),:) + nodes(edges(:,2),:));
facesCenter = 1/4 * ( edgesCenter(faces(:,1),:) + edgesCenter(faces(:,2),:) + ...
                       edgesCenter(faces(:,3),:) + edgesCenter(faces(:,4),:) );
cellsCenter = 1/6 * ( facesCenter(cells(:,1),:) + facesCenter(cells(:,2),:) + ...
                      facesCenter(cells(:,3),:) + facesCenter(cells(:,4),:) + ...
                      facesCenter(cells(:,5),:) + facesCenter(cells(:,6),:) );

% Loop over blocks to make additions
Nblk = size(blkLoc,1);
tol = 0.001; % allow small inaccuracy when locating sheets and lines
for i = 1:Nblk

    xmin = min(blkLoc(i,1:2));
    xmax = max(blkLoc(i,1:2));
    ymin = min(blkLoc(i,3:4));
    ymax = max(blkLoc(i,3:4));
    zmax = max(blkLoc(i,5:6));
    zmin = min(blkLoc(i,5:6));
    [~, xminInd] = min(abs(nodeX-xmin)); % nearest snap to grid
    [~, xmaxInd] = min(abs(nodeX-xmax));
    [~, yminInd] = min(abs(nodeY-ymin));
    [~, ymaxInd] = min(abs(nodeY-ymax));
    [~, zminInd] = min(abs(nodeZ-zmin));
    [~, zmaxInd] = min(abs(nodeZ-zmax));
    
    
    switch objType(i)
        case 3 % volume -> add to cellCon
            
            ind = cellsCenter(:,1) >= nodeX(xminInd) & cellsCenter(:,1) <= nodeX(xmaxInd) & ...
                  cellsCenter(:,2) >= nodeY(yminInd) & cellsCenter(:,2) <= nodeY(ymaxInd) & ...
                  cellsCenter(:,3) <= nodeZ(zmaxInd) & cellsCenter(:,3) >= nodeZ(zminInd);
            cellVal(ind) = blkVal(i);
            
        case 2 % sheet -> add to faceCon
            
            ind = facesCenter(:,1) >= nodeX(xminInd)-tol & facesCenter(:,1) <= nodeX(xmaxInd)+tol & ...
                  facesCenter(:,2) >= nodeY(yminInd)-tol & facesCenter(:,2) <= nodeY(ymaxInd)+tol & ...
                  facesCenter(:,3) <= nodeZ(zmaxInd)+tol & facesCenter(:,3) >= nodeZ(zminInd)-tol;
            faceVal(ind) = blkVal(i);
            
        case 1 % string -> add to edgeCon
            
            ind = edgesCenter(:,1) >= nodeX(xminInd)-tol & edgesCenter(:,1) <= nodeX(xmaxInd)+tol & ...
                  edgesCenter(:,2) >= nodeY(yminInd)-tol & edgesCenter(:,2) <= nodeY(ymaxInd)+tol & ...
                  edgesCenter(:,3) <= nodeZ(zmaxInd)+tol & edgesCenter(:,3) >= nodeZ(zminInd)-tol;
            edgeVal(ind) = blkVal(i);
            
        case 0 % point -> no action
            
    end
    
    
end






end
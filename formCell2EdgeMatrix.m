% Form the mapping matrix that transform cell conductivity model (cellCon in S/m)
% to conductance on edges; mesh-independent.
% function Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes)
% INPUT
%     edges:  a 2-column matrix of node index for the edges; 1st column for
%     starting node and 2nd column for ending node
%     lengths: a vector of the edges' lengths in meter
%     faces: a 4-column matrix of edge index for the faces
%     cells: a 6-column matrix of face index for the cells
%     volumes: a vector of the cells' volume in cubic meter
% OUTPUT
%     Cell2Edge: a Nedges x Ncells matrix that is equivalent to cellVolume/length/length
% NOTE
%     Suppose cellCon is the cell conductivity vector defined on all of the
%     cells. Cc = Cell2Edge * cellCon gets the conductances of the
%     equivalent resistors defined on the edges. The total conductance from
%     a cell is equally distributed to all its edges.
function Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes)

Nedges = size(edges,1);
[Nfaces, Nepf] = size(faces); % # of faces and # of edges per face
[Ncells, Nfpc] = size(cells); % # of cells and # of faces per cell

% form Cell2Edge (matrix to be multiplied by cellCon)
I1 = faces; % row index (which edge)
J1 = repmat((1:Nfaces)',1,Nepf); % column index
S1 = repmat(1/Nepf,Nfaces,Nepf); % % each edge receivers equal portion from face
I2 = cells; % row index (which face)
J2 = repmat((1:Ncells)',1,Nfpc); % column index
S2 = repmat(volumes./Nfpc,1,Nfpc); % equivalent volume per edge - edges equally share the volume from a cell
Cell2Edge = spdiags(1./lengths.^2,0,Nedges,Nedges) * ...
            sparse(I1(:),J1(:),S1(:),Nedges,Nfaces) * ... % distribute volume from face to edge
            sparse(I2(:),J2(:),S2(:),Nfaces,Ncells); % distribute volume from cell to face
        
end


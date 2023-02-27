function Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes)
% A function in the package "RESnet-m" 
% Form the mapping matrix that transforms cell conductivity model (cellCon in S/m)
% to conductance on edges.
%
% function Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes)
% INPUT
%     edges: a 2-column matrix of node index for the edges; 1st column for
%         starting node and 2nd column for ending node
%     lengths: a vector of the edges' lengths in meter
%     faces: a 4-column matrix of edge index for the faces
%     cells: a 6-column matrix of face index for the cells
%     volumes: a vector of the cells' volume in cubic meter
% OUTPUT
%     Cell2Edge: a Nedges x Ncells matrix that acts as volume/4/length/length
% NOTE
%     Suppose cellCon is the cell conductivity vector defined on all of the
%     cells. Cc = Cell2Edge * cellCon gets the conductances of the
%     equivalent resistors defined in the cells.

Nnodes = 8;
[Nedges, ~] = size(edges);
[~, Nepf] = size(faces);
[Ncells, Nfpc] = size(cells);

J = faces( reshape(cells',[],1) , : );
J = reshape(J',[],1);
J = unique( reshape(J,Nepf*Nfpc,[]), 'rows' , 'stable' )';
J = J(:);
I = repmat( (1:Ncells)' , 1 , Nfpc+Nnodes-2 );
I = I(:);
Cell2Edge = spdiags(1./lengths.^2,0,Nedges,Nedges) * sparse(J,I,1,Nedges,Ncells) * spdiags(volumes./4,0,Ncells,Ncells);

end
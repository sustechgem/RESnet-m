function Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas)
% A function in the package "RESnet-m" 
% Form the mapping matrix that transforms face conductivity model (faceCon in S)
% to conductance on edges.
%
% function Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas)
% INPUT
%     edges: a 2-column matrix of node index for the edges; 1st column for
%         starting node and 2nd column for ending node
%     lengths: a vector of the edges' lengths in meter
%     faces: a 4-column matrix of edge index for the faces
%     areas: a vector of the faces' area in square meter
% OUTPUT
%     Face2Edge: a Nedges x Nfaces matrix that acts as area/2/length/length
% NOTE
%     Suppose faceCon is the face conductivity vector defined on all of the
%     faces. Cf = Face2Edge * faceCon gets the conductances of the
%     equivalent resistors defined on the faces.

[Nedges,~] = size(edges); % # of edges
[Nfaces,Nepf] = size(faces); % # of faces, # of edges per face

I = reshape( repmat((1:Nfaces)',1,Nepf)' ,[] ,1 );
J = reshape(faces',[], 1);
Face2Edge = spdiags(1./lengths.^2,0,Nedges,Nedges) * sparse(J,I,1,Nedges,Nfaces) * spdiags(areas./2,0,Nfaces,Nfaces);

end
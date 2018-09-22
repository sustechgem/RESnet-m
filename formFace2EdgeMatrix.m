% Form the mapping matrix that transform face conductivity model (faceCon in S)
% to conductance on edges; mesh-independent.
% function Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas)
% INPUT
%     edges:  a 2-column matrix of node index for the edges; 1st column for
%     starting node and 2nd column for ending node
%     lengths: a vector of the edges' lengths in meter
%     faces: a 4-column matrix of edge index for the faces
%     areas: a vector of the faces' area in square meter
% OUTPUT
%     Face2Edge: a Nedges x Nfaces matrix that is equivalent to faceArea/length/length
% NOTE
%     Suppose faceCon is the face conductivity vector defined on all of the
%     faces. Cf = Face2Edge * faceCon gets the conductances of the
%     equivalent resistors defined on the edges. The total conductance from
%     a face is equally distributed to all its edges.
function Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas)

Nedges = size(edges,1);
[Nfaces, Nepf] = size(faces); % # of faces and # of edges per face

% form Face2Edge (matrix to be multiplied by faceCon)
I = faces; % row index (which edge)
J = repmat((1:Nfaces)',1,Nepf); % column index
S = repmat(areas./Nepf,1,Nepf); % equivalent area per edge - edges equally share the area from a face
Face2Edge = spdiags(1./lengths.^2,0,Nedges,Nedges) * sparse(I(:),J(:),S(:),Nedges,Nfaces); % distribute area from face to edge



end
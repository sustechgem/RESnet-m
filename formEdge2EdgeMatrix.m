function Edge2Edge = formEdge2EdgeMatrix(edges,lengths)
% A function in the package "RESnet-m" 
% Form the mapping matrix that transforms edge conductivity model (edgeCon in S*m)
% to conductance on edges.
%
% function Edge2Edge = formEdge2EdgeMatrix(edges,lengths)
% INPUT
%     edges: a 2-column matrix of node index for the edges; 1st column for
%         starting node and 2nd column for ending node
%     lengths: a vector of the edges' lengths in meter
% OUTPUT
%     Edge2Edge: a Nedges x Nedges matrix that acts as 1/length
% NOTE
%     Suppose edgeCon is the edge conductivity vector defined on all of the
%     edges. Ce = Edge2Edge * edgeCon gets the conductances of the
%     equivalent resistors defined on the edges.

Nedges = size(edges,1);
Edge2Edge = spdiags(1./lengths,0,Nedges,Nedges);

end
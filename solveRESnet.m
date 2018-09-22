% Solve an arbitrary resistor network circuit problem
% function [ potentials, potentialDiffs, currents] = solveRESnet(edges,c,s)
% INPUT
%     edges: a 2-column matrix of node index for the edges that describes 
%     topology of the network; 1st column for starting node and 2nd column 
%     for ending node
%     c: a vector of conductance values on edges
%     s: a vector for the source (current injection amplitude at each node)
% OUTPUT
%     potentials: electric potentials on each node (assume zero potential 
%     on the first node)
%     potentialDiffs: potential drops across each edge
%     currents: current flowing along each edge
function [ potentials, potentialDiffs, currents] = solveRESnet(edges,c,s)

Nnodes = max(max(edges)); % # of nodes
Nedges = size(edges,1); % # of edges

% form potential difference matrix (node to edge), a.k.a. gradient operator
I = kron(1:Nedges,[1; 1]);
J = edges';
S = kron(ones(1,Nedges),[1; -1]);
G = sparse(I(:),J(:),S(:),Nedges,Nnodes); 

C = spdiags(c,0,Nedges,Nedges);
E = sparse(1,1,1,Nnodes,Nnodes,1);
A = G' * C * G + E;

% default solver
potentials = A \ s; 

% compute potential difference (E field) on all edges
potentialDiffs = G * potentials;

% compute current on all edges
currents = C * potentialDiffs;
                         

end
                                     
                                                
   



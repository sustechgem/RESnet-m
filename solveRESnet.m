function [potentials, potentialDiffs, currents] = solveRESnet(edges,C,sources)
% A function in the package "RESnet-m" 
% Solve an arbitrary 3D resistor network circuit problem using the
% potential's formulation and Kirchoff's current law
%
% function [potentials, potentialDiffs, currents] = solveRESnet(edges,C,sources)
% INPUT
%     edges: a 2-column matrix of node index for the edges (branches) that 
%         describes topology of the network; the 1st column for starting node 
%         and the 2nd column for ending node
%     C: a vector of conductance values on edges
%     sources: a vector for the source (current injection amplitude at each node)
% OUTPUT
%     potentials: electric potentials on each node (assume zero potential 
%         at the first node)
%     potentialDiffs: potential drops across each edge (branch)
%     currents: current flowing along each edge (branch)


Nnodes = max(max(edges)); % # of nodes
Nedges = size(edges,1); % # of edges

% Form potential difference matrix (node to edge), a.k.a. gradient operator
I = kron(1:Nedges,[1; 1]);
J = edges';
S = kron(ones(1,Nedges),[1; -1]);
G = sparse(I(:),J(:),S(:),Nedges,Nnodes); 

Cdiag = spdiags(C,0,Nedges,Nedges);
E = sparse(1,1,1,Nnodes,Nnodes,1);
A = G' * Cdiag * G + E;

% Matrix factorization
dA = decomposition(A,'chol');

% Solve for multiple rhs
potentials = dA \ sources; 

% Compute potential difference (E field) on all edges
potentialDiffs = G * potentials;

% Compute current on all edges
currents = Cdiag * potentialDiffs;
                         
end
                                     
                                                
   



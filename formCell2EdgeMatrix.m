function Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes)

Nnodes = 8;
[Nedges,~]    = size(edges);
[~,Nepf] = size(faces);
[Ncells,Nfpc] = size(cells);

J = faces( reshape(cells',[],1) , : );
J = reshape(J',[],1);
J = unique( reshape(J,Nepf*Nfpc,[]), 'rows' , 'stable' )';
J = J(:);
I = repmat( (1:Ncells)' , 1 , Nfpc+Nnodes-2 );
I = I(:);

Cell2Edge = spdiags(1./lengths.^2,0,Nedges,Nedges) * sparse(J,I,1,Nedges,Ncells) * spdiags(volumes./4,0,Ncells,Ncells);

end
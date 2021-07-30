
function Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas)

[Nedges,~] = size(edges);
[Nfaces,Nepf] = size(faces);

I = reshape( repmat((1:Nfaces)',1,Nepf)' ,[] ,1 );
J = reshape(faces',[], 1);

Face2Edge = spdiags(1./lengths.^2,0,Nedges,Nedges) * sparse(J,I,1,Nedges,Nfaces) * spdiags(areas./2,0,Nfaces,Nfaces);

end
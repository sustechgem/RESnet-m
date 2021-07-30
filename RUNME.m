% A quick template for running RESnet
% For demonstration only. Results not numerically meaningful.
% Dikun Yang, 2015 - 2018.

%% Setup the problem

% Create a 3D rectilinear mesh
nodeX = -100:20:100; % node locations in X
nodeY = -100:20:100; % node locations in Y
nodeZ = 0:-20:-500; % node locations in Z

% Define the model; a model includes some blocks that can represent objects like sheets or lines when one or two dimensions vanish.
blkLoc = [ -inf inf -inf inf   0 -inf; % a volumetric object defined by [xmin xmax ymin ymax zmax zmin]
           -20  20  -20  20 -100 -100; % a sheet object defined by [xmin xmax ymin ymax z z]; here is a horizontal sheet
           0    0   0    0   0   -400]; % a line object defined by [x x y y zmax zmin]; here is a vertical line

blkCon = [ 1e-2; % conductive property of the volumetric object (S/m)
           1e1; % conductive property of the sheet object (S)
           5e4]; % conductive property of the line object (S*m)

undefined = 0; % if object value not specified

% Define the current sources [x y z current(A)]; can include multiple lines
source = [0   0   0   1;  % positive current electrode
          100 0   0  -1]; % negative current electrode

%% Form a resistor network

% Get connectivity properties of nodes, edges, faces, cells
[nodes, edges, lengths, faces, areas, cells, volumes] = ...
                                formRectMeshConnectivity(nodeX,nodeY,nodeZ); 

% Get conductive property model vectors (convert the block model description to values on edges, faces and cells)
[edgeCon,faceCon,cellCon] = makeRectMeshModels(nodeX,nodeY,nodeZ,blkLoc,blkCon,undefined);    

% Convert all conductive objects to conductance on edges
Edge2Edge = formEdge2EdgeMatrix(edges,lengths);
Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas);
Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes);
ce = Edge2Edge * edgeCon; % conductance from edges
cf = Face2Edge * faceCon; % conductance from faces
cc = Cell2Edge * cellCon; % conductance from cells
c = ce + cf + cc; % total conductance 

% Source vector (injected current at each node)
s = formSourceNearestNodes(nodes,source); % snap electrodes to the nearest node


%% Solve the network problem

tic
[ potentials, potentialDiffs, currents] = solveRESnet(edges,c,s);                          
toc                            
                            

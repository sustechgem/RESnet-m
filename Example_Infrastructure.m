% Demonstration of complex infrastructure modeling
% Dikun Yang, yangdikun@gmail.com, 2015 - 2023.

% In the following, a dipole-dipole survey is simulated for 
% Model #1: A uniform half-space
% Model #2: Model #1 with two anomalous blocks buried
% Model #3: Model #2 with infrastructure (a steel cased well and a steel-sheet warehouse)
% Model #4: Model #3 with an above-ground pipe connecting the well and the warehouse

%% Setup the 3D mesh

% Create a 3D rectilinear mesh
h = 1; % smallest horizontal grid size 
ratio = 1.4; % expansion rate
nctbc = 15; % number of cells to boundary condition
nodeX = round([ -cumsum(h*ratio.^(nctbc:-1:0)','reverse')-45; (-45:h:45)'; 45+cumsum(h*ratio.^(0:nctbc)')]); % node locations in X
nodeY = round([ -cumsum(h*ratio.^(nctbc:-1:0)','reverse')-12; (-12:h:12)'; 12+cumsum(h*ratio.^(0:nctbc)')]); % node locations in Y

h = 1; % smallest vertical grid size 
ratio = 1.4; % expansion rate
nctbc = 15; % number of cells to boundary condition
nodeZ = [(4:-h:-10)'; -10-round(cumsum(h*ratio.^(0:nctbc)'))]; % node locations in Z (top of mesh at 4 m above the surface)

%% Setup the electric surveys (dipole-dipole)

spacing = 4; % meters between two electrodes
Aloc = -40:spacing:28;
Ntx = length(Aloc); % number of source combos
tx = cell(Ntx,1);
rx = cell(Ntx,1);
for i = 1:Ntx
    tx{i} = [Aloc(i)         0 0  1;  % A electrode
             Aloc(i)+spacing 0 0 -1]; % B electrode
    Mloc = ((Aloc(i)+spacing*2):spacing:36)'; % M electrodes
    Nloc = ((Aloc(i)+spacing*3):spacing:40)'; % N electrodes
    rx{i} = [Mloc Mloc*0 Mloc*0 Nloc Nloc*0 Nloc*0]; % Mx My(=0) Mz(=0) Nx Ny(=0) Nz(=0)
end
Ndata = sum(1:Ntx); % total number of data

%% Define Model #1: Uniform half-space

blkLoc = [ -inf  inf  -inf  inf  inf     0;   % a uniform layer for the air (above surface)
           -inf  inf  -inf  inf    0  -inf];  % a uniform half-space below surface

blkCon = [ 1e-5;   % air conductivity
           0.01];  % earth half-space conductivity

%% Form a resistor network

% Get connectivity properties of nodes, edges, faces, cells
[nodes, edges, lengths, faces, areas, cells, volumes] = ...
                                formRectMeshConnectivity(nodeX,nodeY,nodeZ); 

% Get conductive property model vectors (convert the block-model description to values on edges, faces and cells)
[cellCon,faceCon,edgeCon] = ...
         makeRectMeshModelBlocks(nodeX,nodeY,nodeZ,blkLoc,blkCon,[],[],[]);   

% Convert all conductive objects to conductance on edges
Edge2Edge = formEdge2EdgeMatrix(edges,lengths);
Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas);
Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes);
Ce = Edge2Edge * edgeCon; % conductance from edges
Cf = Face2Edge * faceCon; % conductance from faces
Cc = Cell2Edge * cellCon; % conductance from cells
C = Ce + Cf + Cc; % total conductance 

%% Solve the resistor network problem

% Calculate current sources on the nodes using info in tx
Ntx = length(tx); % number of tx-rx sets
sources = zeros(size(nodes,1),Ntx); % current source intensity at all nodes (for all source configurations)
for i = 1:Ntx
    weights = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,tx{i}(:,1:3)); % weights for the distribution of point current source to the neighboring nodes
    sources(:,i) = weights * tx{i}(:,4); % total current intensities at all the nodes
end

% Obtain potentials at the nodes, potential differences and current along the edges
tic;
[ potentials, potentialDiffs, currents] = solveRESnet(edges,C,sources);                          
toc;                            

% Get simulated data
data = cell(Ntx,1);
for i = 1:Ntx
    Mw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,1:3)); % weights for the interpolation of potential data at the M-electrode location
    Nw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,4:6)); % weights for the interpolation of potential data at the N-electrode location
    data{i} = (Mw' - Nw') * potentials(:,i); % calculate the potential difference data as "M - N"
end

data1 = data; % save to data1

%% Plot Model #1's apparent resistivity

ABMNKVR1 = zeros(Ndata,7); 
count = 0;
for i = 1:Ntx
    A = tx{i}(1,1);
    B = tx{i}(2,1);
    for j = 1:size(rx{i},1)
        M = rx{i}(j,1);
        N = rx{i}(j,4);
        K = 2 * pi / (1/(M-A)-1/(M-B)-1/(N-A)+1/(N-B)); % geometric factor
        count = count + 1;
        ABMNKVR1(count,:) = [A B M N K data1{i}(j) K*data1{i}(j)]; % [A B M N geometric_factor potential_difference apparent_resistivity]
    end
end

x = sum(ABMNKVR1(:,1:4),2)/4; % mid-point of four electrodes
n = (ABMNKVR1(:,3) - ABMNKVR1(:,2))/spacing; % n-spacing as pseudo-depth

figure;
scatter(x,n,200,log10(ABMNKVR1(:,7)),'square','filled');
set(gca,'ydir','reverse');
grid on;
hc = colorbar;
hc.Label.String = 'Apparent resistivity log(Ohm*m)';
xlabel('X (m)');
ylabel('n-spacing');
title('(a) Model #1: Half-space');


%% Define Model #2: Two blocks in half-space

blkLoc = [ -inf  inf  -inf  inf  inf     0;  % a uniform layer for the air (above surface)
           -inf  inf  -inf  inf    0  -inf;  % a uniform half-space below surface
            -30  -20    -4    4   -6   -10;  % an anomalous conductive block
              0    8    -8    2   -2    -6]; % an anomalous resistive block

blkCon = [ 1e-5;   % air conductivity
           0.01;   % earth conductivity
           0.1;    % block conductivity
           0.001]; % block conductivity

%% Form a resistor network

% Get connectivity properties of nodes, edges, faces, cells
[nodes, edges, lengths, faces, areas, cells, volumes] = ...
                                formRectMeshConnectivity(nodeX,nodeY,nodeZ); 

% Get conductive property model vectors (convert the block-model description to values on edges, faces and cells)
[cellCon,faceCon,edgeCon] = ...
         makeRectMeshModelBlocks(nodeX,nodeY,nodeZ,blkLoc,blkCon,[],[],[]);   

% Convert all conductive objects to conductance on edges
Edge2Edge = formEdge2EdgeMatrix(edges,lengths);
Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas);
Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes);
Ce = Edge2Edge * edgeCon; % conductance from edges
Cf = Face2Edge * faceCon; % conductance from faces
Cc = Cell2Edge * cellCon; % conductance from cells
C = Ce + Cf + Cc; % total conductance 

%% Solve the resistor network problem

% Calculate current sources on the nodes using info in tx
Ntx = length(tx); % number of tx-rx sets
sources = zeros(size(nodes,1),Ntx); % current source intensity at all nodes (for all source configurations)
for i = 1:Ntx
    weights = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,tx{i}(:,1:3)); % weights for the distribution of point current source to the neighboring nodes
    sources(:,i) = weights * tx{i}(:,4); % total current intensities at all the nodes
end

% Obtain potentials at the nodes, potential differences and current along the edges
tic;
[ potentials, potentialDiffs, currents] = solveRESnet(edges,C,sources);                          
toc;                            

% Get simulated data
data = cell(Ntx,1);
for i = 1:Ntx
    Mw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,1:3)); % weights for the interpolation of potential data at the M-electrode location
    Nw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,4:6)); % weights for the interpolation of potential data at the N-electrode location
    data{i} = (Mw' - Nw') * potentials(:,i); % calculate the potential difference data as "M - N"
end

data2 = data; % save to data2

%% Plot Model #2's apparent resistivity

ABMNKVR2 = zeros(Ndata,7); 
count = 0;
for i = 1:Ntx
    A = tx{i}(1,1);
    B = tx{i}(2,1);
    for j = 1:size(rx{i},1)
        M = rx{i}(j,1);
        N = rx{i}(j,4);
        K = 2 * pi / (1/(M-A)-1/(M-B)-1/(N-A)+1/(N-B)); % geometric factor
        count = count + 1;
        ABMNKVR2(count,:) = [A B M N K data2{i}(j) K*data2{i}(j)]; % [A B M N geometric_factor potential_difference apparent_resistivity]
    end
end

x = sum(ABMNKVR2(:,1:4),2)/4; % mid-point of four electrodes
n = (ABMNKVR2(:,3) - ABMNKVR2(:,2))/spacing; % n-spacing as pseudo-depth

figure;
scatter(x,n,200,log10(ABMNKVR2(:,7)),'square','filled');
set(gca,'ydir','reverse');
grid on;
hc = colorbar;
hc.Label.String = 'Apparent resistivity log(Ohm*m)';
xlabel('X (m)');
ylabel('n-spacing');
title('(b) Model #2: Two blocks');

%% Define Model #3: Infrastructure

blkLoc = [ -inf  inf  -inf  inf  inf     0;  % a uniform layer for the air (above surface)
           -inf  inf  -inf  inf    0  -inf;  % a uniform half-space below surface
            -30  -20    -4    4   -6   -10;  % an anomalous conductive block
              0    8    -8    2   -2    -6;  % an anomalous resistive block
            -18  -18    -6   -6    0  -100;  % a steel cased well
              6   16     4    4    6    -2;  % southern wall of steel-sheet warehouse
              6   16    10   10    6    -2;  % northern wall of steel-sheet warehouse
              6    6     4   10    6    -2;  % western wall of steel-sheet warehouse
             16   16     4   10    6    -2;  % eastern wall of steel-sheet warehouse
              6   16     4   10    6     6]; % roof

blkCon = [ 1e-5;   % air conductivity
           0.01;   % earth conductivity
           0.1;    % block conductivity
           0.001;  % block conductivity
           pi*(0.05^2-0.04^2) * 5e6 ;  % steel casing's edgeCon (cross-sectional area times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6]; % steel sheet's faceCon (thickness times intrinsic conductivity)

 
%% Form a resistor network

% Get connectivity properties of nodes, edges, faces, cells
[nodes, edges, lengths, faces, areas, cells, volumes] = ...
                                formRectMeshConnectivity(nodeX,nodeY,nodeZ); 

% Get conductive property model vectors (convert the block-model description to values on edges, faces and cells)
[cellCon,faceCon,edgeCon] = ...
         makeRectMeshModelBlocks(nodeX,nodeY,nodeZ,blkLoc,blkCon,[],[],[]);   

% Convert all conductive objects to conductance on edges
Edge2Edge = formEdge2EdgeMatrix(edges,lengths);
Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas);
Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes);
Ce = Edge2Edge * edgeCon; % conductance from edges
Cf = Face2Edge * faceCon; % conductance from faces
Cc = Cell2Edge * cellCon; % conductance from cells
C = Ce + Cf + Cc; % total conductance 

%% Solve the resistor network problem

% Calculate current sources on the nodes using info in tx
Ntx = length(tx); % number of tx-rx sets
sources = zeros(size(nodes,1),Ntx); % current source intensity at all nodes (for all source configurations)
for i = 1:Ntx
    weights = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,tx{i}(:,1:3)); % weights for the distribution of point current source to the neighboring nodes
    sources(:,i) = weights * tx{i}(:,4); % total current intensities at all the nodes
end

% Obtain potentials at the nodes, potential differences and current along the edges
tic;
[ potentials, potentialDiffs, currents] = solveRESnet(edges,C,sources);                          
toc;                            

% Get simulated data
data = cell(Ntx,1);
for i = 1:Ntx
    Mw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,1:3)); % weights for the interpolation of potential data at the M-electrode location
    Nw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,4:6)); % weights for the interpolation of potential data at the N-electrode location
    data{i} = (Mw' - Nw') * potentials(:,i); % calculate the potential difference data as "M - N"
end

data3 = data; % save to data3

%% Plot Model #3's apparent resistivity

ABMNKVR3 = zeros(Ndata,7); 
count = 0;
for i = 1:Ntx
    A = tx{i}(1,1);
    B = tx{i}(2,1);
    for j = 1:size(rx{i},1)
        M = rx{i}(j,1);
        N = rx{i}(j,4);
        K = 2 * pi / (1/(M-A)-1/(M-B)-1/(N-A)+1/(N-B)); % geometric factor
        count = count + 1;
        ABMNKVR3(count,:) = [A B M N K data3{i}(j) K*data3{i}(j)]; % [A B M N geometric_factor potential_difference apparent_resistivity]
    end
end

x = sum(ABMNKVR3(:,1:4),2)/4; % mid-point of four electrodes
n = (ABMNKVR3(:,3) - ABMNKVR3(:,2))/spacing; % n-spacing as pseudo-depth

figure;
scatter(x,n,200,log10(ABMNKVR3(:,7)),'square','filled');
set(gca,'ydir','reverse');
grid on;
hc = colorbar;
hc.Label.String = 'Apparent resistivity log(Ohm*m)';
xlabel('X (m)');
ylabel('n-spacing');
title('(c) Model #3: Infrastructure');


%% Define Model #4: Above-ground pipe

blkLoc = [ -inf  inf  -inf  inf  inf     0;  % a uniform layer for the air (above surface)
           -inf  inf  -inf  inf    0  -inf;  % a uniform half-space below surface
            -30  -20    -4    4   -6   -10;  % an anomalous conductive block
              0    8    -8    2   -2    -6;  % an anomalous resistive block
            -18  -18    -6   -6    4  -100;  % a steel cased well
              6   16     4    4    6    -2;  % southern wall of steel-sheet warehouse
              6   16    10   10    6    -2;  % northern wall of steel-sheet warehouse
              6    6     4   10    6    -2;  % western wall of steel-sheet warehouse
             16   16     4   10    6    -2;  % eastern wall of steel-sheet warehouse
              6   16     4   10    6     6]; % roof

blkCon = [ 1e-5;   % air conductivity
           0.01;   % earth conductivity
           0.1;    % block conductivity
           0.001;  % block conductivity
           pi*(0.05^2-0.04^2) * 5e6 ;  % steel casing's edgeCon (cross-sectional area times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6;  % steel sheet's faceCon (thickness times intrinsic conductivity)
           0.005 * 5e6]; % steel sheet's faceCon (thickness times intrinsic conductivity)
 
%% Form a resistor network

% Get connectivity properties of nodes, edges, faces, cells
[nodes, edges, lengths, faces, areas, cells, volumes] = ...
                                formRectMeshConnectivity(nodeX,nodeY,nodeZ); 

% Get conductive property model vectors (convert the block-model description to values on edges, faces and cells)
[cellCon,faceCon,edgeCon] = ...
         makeRectMeshModelBlocks(nodeX,nodeY,nodeZ,blkLoc,blkCon,[],[],[]);   

% Convert all conductive objects to conductance on edges
Edge2Edge = formEdge2EdgeMatrix(edges,lengths);
Face2Edge = formFace2EdgeMatrix(edges,lengths,faces,areas);
Cell2Edge = formCell2EdgeMatrix(edges,lengths,faces,cells,volumes);
Ce = Edge2Edge * edgeCon; % conductance from edges
Cf = Face2Edge * faceCon; % conductance from faces
Cc = Cell2Edge * cellCon; % conductance from cells
C = Ce + Cf + Cc; % total conductance 

% Add above-ground pipe
pipeStart = [-18, -6, 4]; % where the pipe starts
pipeEnd = [6, 4, 4];      % where the pipe ends
pipeLength = norm(pipeStart-pipeEnd); % length of the pipe
pipeStartNode = find(nodes(:,1)==pipeStart(1) & nodes(:,2)==pipeStart(2) & nodes(:,3)==pipeStart(3)); % search for the starting node
pipeEndNode = find(nodes(:,1)==pipeEnd(1) & nodes(:,2)==pipeEnd(2) & nodes(:,3)==pipeEnd(3));  % search for the ending node
edges = [edges; pipeStartNode pipeEndNode]; % append an additional edge representing the pipe 
C = [C; pi*(0.05^2-0.04^2)*5e6/pipeLength]; % append an additional conductance representing the pipe


%% Solve the resistor network problem

% Calculate current sources on the nodes using info in tx
Ntx = length(tx); % number of tx-rx sets
sources = zeros(size(nodes,1),Ntx); % current source intensity at all nodes (for all source configurations)
for i = 1:Ntx
    weights = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,tx{i}(:,1:3)); % weights for the distribution of point current source to the neighboring nodes
    sources(:,i) = weights * tx{i}(:,4); % total current intensities at all the nodes
end

% Obtain potentials at the nodes, potential differences and current along the edges
tic;
[ potentials, potentialDiffs, currents] = solveRESnet(edges,C,sources);                          
toc;                            

% Get simulated data
data = cell(Ntx,1);
for i = 1:Ntx
    Mw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,1:3)); % weights for the interpolation of potential data at the M-electrode location
    Nw = calcTrilinearInterpWeights(nodeX,nodeY,nodeZ,rx{i}(:,4:6)); % weights for the interpolation of potential data at the N-electrode location
    data{i} = (Mw' - Nw') * potentials(:,i); % calculate the potential difference data as "M - N"
end

data4 = data; % save to data4

%% Plot Model #4's apparent resistivity

ABMNKVR4 = zeros(Ndata,7); 
count = 0;
for i = 1:Ntx
    A = tx{i}(1,1);
    B = tx{i}(2,1);
    for j = 1:size(rx{i},1)
        M = rx{i}(j,1);
        N = rx{i}(j,4);
        K = 2 * pi / (1/(M-A)-1/(M-B)-1/(N-A)+1/(N-B)); % geometric factor
        count = count + 1;
        ABMNKVR4(count,:) = [A B M N K data4{i}(j) K*data4{i}(j)]; % [A B M N geometric_factor potential_difference apparent_resistivity]
    end
end

x = sum(ABMNKVR4(:,1:4),2)/4; % mid-point of four electrodes
n = (ABMNKVR4(:,3) - ABMNKVR4(:,2))/spacing; % n-spacing as pseudo-depth

figure;
scatter(x,n,200,log10(ABMNKVR4(:,7)),'square','filled');
set(gca,'ydir','reverse');
grid on;
hc = colorbar;
hc.Label.String = 'Apparent resistivity log(Ohm*m)';
xlabel('X (m)');
ylabel('n-spacing');
title('(d) Model #4: Above-ground pipe');





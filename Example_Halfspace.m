% Testing the solution of half-space
% Dikun Yang, yangdikun@gmail.com, 2015 - 2023.

% In the following, a pole-dipole array is simulated for a uniform
% half-space and compared against the analytic solution.

%% Setup the 3D mesh

% Create a 3D rectilinear mesh
h = 2;
ratio = 1.14;
nctbc = 30;
nodeX = round([ -cumsum(h*ratio.^(nctbc:-1:0)','reverse'); 0; cumsum(h*ratio.^(0:nctbc)')]); % node locations in X
nodeY = round([ -cumsum(h*ratio.^(nctbc:-1:0)','reverse'); 0; cumsum(h*ratio.^(0:nctbc)')]); % node locations in Y
nodeZ = round([0; -cumsum(h*ratio.^(0:nctbc)')]); % node locations in Z

%% Setup the geo-electrical model

% Define the model using combination of blocks
blkLoc = [ -inf  inf  -inf  inf    0 -inf];  % a uniform half-space
blkCon = [ 1e-2]; % conductive property of a volumetric object (S/m)

%% Setup the electric surveys (pole-dipole)

% Define the current sources in the format of [x y z current(Ampere)]
tx = {[    0    0    0    1;   % A electrode
        -inf    0    0   -1]}; % B electrode

% Define the receiver electrodes in the format of [Mx My Mz Nx Ny Nz]
rx = {[   10    0    0   20    0    0;  % nine M-N pairs for the source
          20    0    0   30    0    0;
          30    0    0   40    0    0;
          40    0    0   50    0    0;
          50    0    0   60    0    0;
          60    0    0   70    0    0;
          70    0    0   80    0    0;
          80    0    0   90    0    0;
          90    0    0  100    0    0]};

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
sources = zeros(size(nodes,1),Ntx);
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

%% Compare against analytic solutions

Aloc = [0, 0, 0]; % location of A electrode
rAM = rx{1}(:,1) - Aloc(1); % A-M distance
rAN = rx{1}(:,4) - Aloc(1); % A-N distance
rho = 100; % half-space resistivity
I = 1;
dV = rho * I / 2 / pi * (1./rAM - 1./rAN); % potential differences (analytic solution)
X = 0.5 * (rx{1}(:,1) + rx{1}(:,4)); % centers of M-N (x-coordinate)

figure;
subplot(2,1,1);
semilogy(X,data{1},'.-');
hold on;
plot(X,dV,'.-');
legend({'RESnet','Analytic'});
grid on;
xlabel('Tx-Rx offset (m)');
ylabel('Potential difference (V)');
title('(a) Numerical and analytic solutions');

subplot(2,1,2);
plot(X,(data{1}-dV)./dV,'k.-');
grid on;
xlabel('Tx-Rx offset (m)');
ylabel('Relative error');
title('(b) Numerical errors');

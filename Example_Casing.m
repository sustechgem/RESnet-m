% Demonstration of steel casing modeling
% Dikun Yang, yangdikun@gmail.com, 2015 - 2023.

% This script reproduces the top row panels (a, b, c, d) in Figure 5 in
% Heagy, L.J. and Oldenburg, D.W., 2019. Direct current resistivity with 
% steel-cased wells. Geophysical Journal International, 219(1), pp.1-26. 
% doi:10.1093/gji/ggz281

%% Setup the geo-electrical model

% Create a 3D rectilinear mesh
h = 50; % smallest horizontal grid size 
ratio = 1.4; % expansion rate
nctbc = 11; % number of cells to boundary condition
nodeX = round([ -cumsum(h*ratio.^(nctbc:-1:0)','reverse')-1300; (-1300:h:1300)'; 1300+cumsum(h*ratio.^(0:nctbc)')]); % node locations in X
nodeY = round([ -cumsum(h*ratio.^(nctbc:-1:0)','reverse')-1300; (-1300:h:1300)'; 1300+cumsum(h*ratio.^(0:nctbc)')]); % node locations in Y

h = 50; % smallest vertical grid size 
ratio = 1.4; % expansion rate
nctbc = 12; % number of cells to boundary condition
nodeZ = [(0:-h:-1000)'; -1000-round(cumsum(h*ratio.^(0:nctbc)'))]; % node locations in Z

% Define the model using combination of blocks
% A model includes some blocks that can represent objects like sheets or lines when one or two dimensions vanish.
blkLoc = [ -inf  inf  -inf  inf    0  -inf;  % a uniform half-space
              0    0     0    0    0 -1000]; % a steel casing


blkCon = [ 0.1; % conductive property of the half-space earth (S/m)
           pi*(0.05^2-0.04^2) * 5e6 ];  % conductive property of the casing (S*m) = cross-sectional area (m^2) * steel's conductivity (S/m)

%% Setup the electric surveys

% Define the current sources in the format of [x y z current(Ampere)]
tx = {[    0    0    0    1;   % first set of source (two electrodes)
       -2000    0    0   -1];  %   return electrode 2000 m away
      [    0    0    0    1;   % second set of source (two electrodes)
        -750    0    0   -1];  %   return electrode 750 m away
      [    0    0    0    1;   % third set of source (two electrodes)
        -500    0    0   -1];  %   return electrode 500 m away
      [    0    0    0    1;   % forth set of source (two electrodes)
        -250    0    0   -1]}; %   return electrode 250 m away

% Define the receiver electrodes in the format of [Mx My Mz Nx Ny Nz]
spacing = 20; % M-N distance
datagridx = -1250:spacing:1250; % data grid in X
Ndatagridx = length(datagridx); % number of data grid in X
datagridy = -1250:spacing:1250; % data grid in Y
Ndatagridy = length(datagridy); % number of data grid in Y
N = Ndatagridx * Ndatagridy; % total number of receivers
[datax, datay] = meshgrid(datagridx,datagridy);
rx = [datax(:)-spacing/2  datay(:)            zeros(N,1)  datax(:)+spacing/2  datay(:)            zeros(N,1);  % measuring Ex
      datax(:)            datay(:)-spacing/2  zeros(N,1)  datax(:)            datay(:)+spacing/2  zeros(N,1)]; % measuring Ey

rx = {rx;  % same receiver locations for four different source configurations
      rx; 
      rx; 
      rx};

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

%% Plot the results

% The top row in Figure 5 (Heagy and Oldenburg, 2019)
figure;
E1 = data{1}/spacing; % first set of source
Ex1 = reshape(E1(1:N),Ndatagridy,Ndatagridx);
Ey1 = reshape(E1(N+1:2*N),Ndatagridy,Ndatagridx);
Etotal1 = sqrt(Ex1.^2+Ey1.^2);
imagesc(datagridx,datagridy,log10(Etotal1));
set(gca,'ydir','reverse');
hc1 = colorbar;
hc1.Label.String = 'primary electric field log(V/m)';
caxis([-7 -4]);
hold on;
contour(datax,datay,log10(Etotal1),'black','ShowText','on');
grid on;
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
title('(a) Return electrode offset = 2000 m');

figure;
E2 = data{2}/spacing; % second set of source
Ex2 = reshape(E2(1:N),Ndatagridy,Ndatagridx);
Ey2 = reshape(E2(N+1:2*N),Ndatagridy,Ndatagridx);
Etotal2 = sqrt(Ex2.^2+Ey2.^2);
imagesc(datagridx,datagridy,log10(Etotal2));
set(gca,'ydir','reverse');
hc2 = colorbar;
hc2.Label.String = 'primary electric field log(V/m)';
caxis([-7 -4]);
hold on;
contour(datax,datay,log10(Etotal2),'black','ShowText','on');
grid on;
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
title('(b) Return electrode offset = 750 m');

figure;
E3 = data{3}/spacing; % third set of source
Ex3 = reshape(E3(1:N),Ndatagridy,Ndatagridx);
Ey3 = reshape(E3(N+1:2*N),Ndatagridy,Ndatagridx);
Etotal3 = sqrt(Ex3.^2+Ey3.^2);
imagesc(datagridx,datagridy,log10(Etotal3));
set(gca,'ydir','reverse');
hc3 = colorbar;
hc3.Label.String = 'primary electric field log(V/m)';
caxis([-7 -4]);
hold on;
contour(datax,datay,log10(Etotal3),'black','ShowText','on');
grid on;
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
title('(c) Return electrode offset = 500 m');

figure;
E4 = data{4}/spacing; % forth set of source
Ex4 = reshape(E4(1:N),Ndatagridy,Ndatagridx);
Ey4= reshape(E4(N+1:2*N),Ndatagridy,Ndatagridx);
Etotal4 = sqrt(Ex4.^2+Ey4.^2);
imagesc(datagridx,datagridy,log10(Etotal4));
set(gca,'ydir','reverse');
hc4 = colorbar;
hc4.Label.String = 'primary electric field log(V/m)';
caxis([-7 -4]);
hold on;
contour(datax,datay,log10(Etotal4),'black','ShowText','on');
grid on;
axis equal;
xlabel('X (m)');
ylabel('Y (m)');
title('(d) Return electrode offset = 250 m');




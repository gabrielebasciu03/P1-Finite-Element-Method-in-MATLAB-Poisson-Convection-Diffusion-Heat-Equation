function [Nelem, Nodes, triang, coord, NDir, DirNod, DirVal, ureal, force,NeuNod] = input_data(meshdir)
% INPUT_DATA  Read mesh files and problem data for the FEM Poisson solver.
%
% Files expected in meshdir/:
%   triang.dat  - element/node connectivity (Nel x 4, last col = zone id)
%   xy.dat      - node coordinates (Nodes x 2)
%   dirnod.dat  - indices of Dirichlet boundary nodes
%   neunod.dat  - indices of Neumann boundary nodes

meshtriang = fullfile(meshdir,'triang.dat');
meshcoord  = fullfile(meshdir,'xy.dat');
meshDirNod = fullfile(meshdir,'dirnod.dat');
meshNeuNod = fullfile(meshdir,'neunod.dat');

triang = load(meshtriang);
coord  = load(meshcoord);


triang = triang(:,1:3);          % drop the zone-id column


Nelem = size(triang,1);
Nodes = size(coord,1);


DirNod = load(meshDirNod);
DirNod = DirNod(:);        % forza sempre vettore colonna
NDir = numel(DirNod);      % conta gli elementi, non le righe
DirVal = zeros(NDir, 1);
NeuNod = load(meshNeuNod);
NeuNod = NeuNod(:);
% Exact solution and forcing term (Section 1 of the project).
% NOTE: f = -Laplacian(ureal); the sign here is the one that is
% consistent with the given exact solution u=sin(pi x)cos(pi y/2).
ureal = @(x,y) sin(pi.*x) .* cos(pi.*y/2);
force = @(x,y)  (5/4) * pi^2 * ureal(x,y);
end

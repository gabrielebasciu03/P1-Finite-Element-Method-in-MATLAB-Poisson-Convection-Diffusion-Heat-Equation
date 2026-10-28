%%% Finite element code for 2D piecewise linear Galekrin 
%%% -div(Diff grad u)(x)=f(x) with Dirichlet BC 
clear all; close all; clc;

% Define the 2D mesh and equation BC's and coefficients
meshdir='mesh2';
[Nelem,Nodes,triang,coord,NDir,DirNod,ureal,force]=input_data(meshdir);

% calculate element area and elemental coefficients of basis functions
% (b,c)
[Bloc,Cloc,Area]=localBasis(Nelem,triang,coord);

% build stiffness matrix (without BCs)
stiffMat=stiffBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,Diff);

% build convection matrix (without BCs)
convMat=convBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,beta);

% build stabilization matrix (without BCs)
stabMat=stabBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,delta);

 sysMat=stiffMat+convMat+stabMat;
% impose BCs
[sysMat,rhs]=imposeBC(Nodes,NDir,DirNod,DirVal,sysMat);

%uh=sysMat\rhs;
ur=ureal(coord(:,1),coord(:,2));
trisurf(triang,coord(:,1),coord(:,2),ur);
 plot_tri(coord,triang,1,'node_index_l')
 plot_tri(coord,triang,1,'edges')
 plot_tri(coord,triang,1,'face_index')
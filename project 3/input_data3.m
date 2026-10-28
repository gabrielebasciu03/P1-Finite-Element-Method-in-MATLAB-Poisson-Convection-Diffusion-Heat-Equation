function [Nelem,Nodes,triang,coord,NDir,DirNod] = input_data3(meshdir)
% Load mesh data and build the FULL Dirichlet boundary (dirnod U neunod),
% since Project 3 requires u=0 on the whole boundary d(Omega).

    load( strcat(meshdir,'/triang.dat') );  % -> triang
    load( strcat(meshdir,'/xy.dat')     );  % -> xy

    triang = triang(:,1:3);
    coord  = xy(:,1:2);
    Nelem  = size(triang,1);
    Nodes  = size(coord,1);

    dirnod_file = strcat(meshdir,'/dirnod.dat');
    dirtot_file = strcat(meshdir,'/dirtot.dat');
    neunod_file = strcat(meshdir,'/neunod.dat');

    if exist(dirnod_file,'file')
        load(dirnod_file);          % -> dirnod
        bnd = dirnod(:);
    elseif exist(dirtot_file,'file')
        load(dirtot_file);          % -> dirtot
        bnd = dirtot(:);
    else
        bnd = [];
    end

    if exist(neunod_file,'file')
        load(neunod_file);          % -> neunod
        bnd = [bnd; neunod(:)];
    end

    DirNod = unique(bnd);
    NDir   = length(DirNod);
end
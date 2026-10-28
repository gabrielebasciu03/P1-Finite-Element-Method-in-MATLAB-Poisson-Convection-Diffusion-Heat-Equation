function [Nelem, Nodes, triang, coord, NDir, DirNod, DirVal, NeuNod] = input_data2(meshdir)
% INPUT_DATA2  Reads the mesh and builds the boundary conditions for
% the convection-diffusion problem of Project 2.
%
% Project 2 instead defines the boundary conditions in a
% purely GEOMETRIC way:
%
%         Gamma1 (u = 1):  [x = 0, 0<=y<=1]  UNION  [0<=x<=0.3, y = 0]
%         Gamma2 (u = 0):  [x = 1, 0<=y<=1]  UNION  [0<=x<=1, y = 1]
%                                            UNION  [0.3<=x<=1, y = 0]
%
% so here we do NOT read dirnod.dat/neunod.dat from disk: we find
% the boundary nodes ourselves (based on their coordinates) and
% classify them into Gamma1/Gamma2. ALL of the domain boundary is of
% Dirichlet type in this problem, so NeuNod is returned
% empty -- rhsBuild2/rhsBuild will not add any flux
% contribution.
% There is no exact solution nor forcing term: the equation
% is homogeneous (f = 0). The forcing is therefore set to zero in the
% main script.
%
% Expected files in meshdir/:
%   triang.dat  - element/node connectivity (Nel x 4, last column = zone)
%   xy.dat      - node coordinates (Nodes x 2)

    meshtriang = fullfile(meshdir, 'triang.dat');
    meshcoord  = fullfile(meshdir, 'xy.dat');

    triang = load(meshtriang);
    coord  = load(meshcoord);

    triang = triang(:, 1:3);   % drop della colonna zone-id

    Nelem = size(triang, 1);
    Nodes = size(coord, 1);

    x = coord(:, 1);
    y = coord(:, 2);

% --- IMPORTANT ---
    % The Project 2 text explicitly requires Omega = [0,1]x[0,1].
% The available meshes, however, come from Project 1
% and are defined on [-1,1]x[-1,1] (verified: min/max of the
% coordinates return -1 and 1).
%
% To solve the problem exactly as specified by the text,
% instead of adapting the Gamma1/Gamma2 formulas to the "wrong" domain,
% we rescale here the mesh coordinates onto the unit square
% [0,1]x[0,1] using an affine transformation (translate + scale). This
% operation does NOT alter the mesh topology (which nodes
% belong to which triangles remains identical, see "triang"),
% it only changes the node coordinates -- exactly as if the mesh
% had been generated directly on the correct domain.   
    xmin = min(x); xmax = max(x);
    ymin = min(y); ymax = max(y);
    x = (x - xmin) / (xmax - xmin);
    y = (y - ymin) / (ymax - ymin);
    coord = [x, y];

    tol = 1e-8;

    % A node is a boundary node if it lies on one of the 4 sides of the square [0,1]x[0,1]
    onBoundary = (abs(x - 0) < tol) | (abs(x - 1) < tol) | ...
                 (abs(y - 0) < tol) | (abs(y - 1) < tol);

    DirNod = find(onBoundary);
    NDir   = numel(DirNod);
    DirVal = zeros(NDir, 1);

    for k = 1:NDir
        idx = DirNod(k);
        xi = x(idx);
        yi = y(idx);

        if abs(xi - 0) < tol
            % Gamma1: side x = 0 (also includes the corner (0,1), which is
            % a point of discontinuity of the boundary data between
           
            DirVal(k) = 1;
        elseif (abs(yi - 0) < tol) && (xi <= 0.3 + tol)
            % Gamma1: edge y = 0, 0 <= x <= 0.3
            DirVal(k) = 1;
        else
            % Gamma2: x = 1, or y = 1, or y = 0 with x > 0.3
            DirVal(k) = 0;
        end
    end

    NeuNod = [];  % no Neumann boundary in this problem

    % Controllo diagnostico: avvisa se pochi nodi risultano di bordo
    % rispetto a quanto atteso (sintomo di un dominio non rilevato bene)
    if NDir < 3
        warning('input_data2:FewNodesBoundary', ...
            'Found only %d boundary nodes out of %d total nodes: check that xy.dat actually contains the domain boundary.', NDir, Nodes);
    end
end
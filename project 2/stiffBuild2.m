function stiffMat = stiffBuild2(Nelem, Nodes, triang, Bloc, Cloc, Area, D)
% STIFFBUILD  Assembles the global stiffness (diffusion) matrix.
%
%
%   If D is not passed, it defaults to 1 and the function behaves
%   exactly as in Project 1 (useful for the comparison required by the
%   text: "first try the code with vx = vy = 0 and compare its solution
%   with the solution obtained with the available Poisson solver").

    if nargin < 7 || isempty(D)
        D = 1;
    end
    if isscalar(D)
        D = D * ones(Nelem, 1);
    end

    stiffMat = sparse(Nodes, Nodes);

    for el = 1:Nelem
        % Local stiffness matrix 3x3
        Sloc = zeros(3, 3);
        for i = 1:3
            for j = 1:3
                Sloc(i,j) = (D(el) / (4 * Area(el))) * (Bloc(el,i)*Bloc(el,j) + Cloc(el,i)*Cloc(el,j));
            end
        end
        if Area(el) <= 0
            error('Not valid Area in element %d. Check the meshes!', el);
        end

        % Assembling
        for i = 1:3
            iglob = triang(el, i);
            for j = 1:3
                jglob = triang(el, j);
                stiffMat(iglob, jglob) = stiffMat(iglob, jglob) + Sloc(i, j);
            end
        end
    end
end

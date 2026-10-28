function [stiffMat, rhs] = imposeBC2(NDir, DirNod, DirVal, stiffMat, rhs, eps_pen)
% IMPOSEBC2  Imposition of Dirichlet boundary conditions via the penalty
% method with penalty value eps = 1e-10.

    if nargin < 6 || isempty(eps_pen)
        eps_pen = 1e-10;
    end

    for i = 1:NDir
        idx = DirNod(i);

        % Add a huge value (1/eps) to the diagonal of the Dirichlet node
        stiffMat(idx, idx) = stiffMat(idx, idx) + 1/eps_pen;

        % Update the rhs consistently with the same penalty
        rhs(idx) = rhs(idx) + DirVal(i) / eps_pen;
    end
end

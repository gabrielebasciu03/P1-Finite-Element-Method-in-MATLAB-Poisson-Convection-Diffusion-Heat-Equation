function [stiffMat, rhs] = imposeBC(NDir, DirNod, DirVal, stiffMat, rhs)
% IMPOSEBC  Imposition of Dirichlet boundary conditions via the penalty
% method with penalty value eps = 1e-15.
%
% Idea: the equation for the Dirichlet node is replaced with
%   (1/eps) * u_idx = DirVal(i) / eps
% which, since eps is very small, numerically "forces" u_idx = DirVal(i)
% while keeping the matrix symmetric.
    eps_pen = 1e-15;

    for i = 1:NDir
        idx = DirNod(i);

        % Add a huge value (1/eps) to the diagonal of the Dirichlet node
        stiffMat(idx, idx) = stiffMat(idx, idx) + 1/eps_pen;

        % Update the rhs consistently with the same penalty
        rhs(idx) = rhs(idx) + DirVal(i) / eps_pen;
    end
end
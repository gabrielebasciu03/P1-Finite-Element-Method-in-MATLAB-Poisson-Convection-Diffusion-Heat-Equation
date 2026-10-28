function transMat = transportBuild(Nelem, Nodes, triang, Bloc, Cloc, Area, betax, betay)
% TRANSPORTBUILD  Assembles the global transport (convection) matrix B.
%
%   b_ij = Integral_Omega (betax * dphi_j/dx + betay * dphi_j/dy) * phi_i dx
%
% Note: unlike the stiffness matrix A, B is NOT symmetric
% (B ~= B^T), as explicitly observed in the project notes.
%
% beta = (betax, betay) is the velocity field, here assumed constant over
% the whole domain (as required by Project 2, beta = (1,3)^T).
%
% Derivation of the local formula:
%   The gradient of phi_j is constant on the triangle T:
%       dphi_j/dx = Bloc(el,j) / (2*Area(el))
%       dphi_j/dy = Cloc(el,j) / (2*Area(el))
%   For P1 elements, Integral_T phi_i dx = Area(el)/3 holds for every i.
%   Therefore:
%       b_ij^T = (betax*Bloc(el,j) + betay*Cloc(el,j)) / (2*Area(el)) * Area(el)/3
%              = (betax*Bloc(el,j) + betay*Cloc(el,j)) / 6
%   which does NOT depend on i: column j of the local matrix has the same
%   value repeated over the 3 rows (a direct consequence of Int_T phi_i = A/3
%   being equal for i=1,2,3).

    transMat = sparse(Nodes, Nodes);

    for el = 1:Nelem
        Tloc = zeros(3, 3);
        for j = 1:3
            val = (betax * Bloc(el, j) + betay * Cloc(el, j)) / 6;
            Tloc(:, j) = val; % same value fot the three rows (i = 1,2,3)
        end

        % Assembling
        for i = 1:3
            iglob = triang(el, i);
            for j = 1:3
                jglob = triang(el, j);
                transMat(iglob, jglob) = transMat(iglob, jglob) + Tloc(i, j);
            end
        end
    end
end

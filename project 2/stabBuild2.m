function stabMat = stabBuild2(Nelem, Nodes, triang, Bloc, Cloc, Area, betax, betay, delta)
% STABBUILD2  Assembles the global Streamline-Upwind Diffusion (SUD)
% stabilization matrix S.
%   s_ij^T = delta * Integral_T (beta.grad(phi_i)) * (beta.grad(phi_j)) dx
%
% with beta = (betax, betay) constant and delta an empirical stabilization
% coefficient
%
% Since the gradient of phi_i and phi_j is constant on the triangle T:
%   beta.grad(phi_i) = (betax*Bloc(el,i) + betay*Cloc(el,i)) / (2*Area(el))
% and the integral of a constant over T equals Area(el) times the constant:
%   s_ij^T = delta * Area(el) *
%            [(betax*Bloc(el,i)+betay*Cloc(el,i)) / (2*Area(el))] *
%            [(betax*Bloc(el,j)+betay*Cloc(el,j)) / (2*Area(el))]
%          = delta / (4*Area(el)) *
%            (betax*Bloc(el,i)+betay*Cloc(el,i)) *
%            (betax*Bloc(el,j)+betay*Cloc(el,j))
%
% Note: S is symmetric and positive semi-definite (it is a "stiffness-type"
% matrix in the direction of beta), and it is added to the system
% matrix A+B exactly as an additional diffusive contribution.
%
% delta can be passed either as a constant scalar
% or as an Nelem x 1 vector if one wants to use a
% "physical" variant delta_T = delta0 * hT / (2*|beta|), with hT the length
% of the longest side of triangle T.

    if isscalar(delta)
        delta = delta * ones(Nelem, 1);
    end

    betaNorm2 = betax^2 + betay^2;
    if betaNorm2 == 0
        stabMat = sparse(Nodes, Nodes);
        return
    end

    stabMat = sparse(Nodes, Nodes);

    for el = 1:Nelem
        Sloc = zeros(3, 3);
        gradDotBeta = zeros(3, 1);
        for i = 1:3
            gradDotBeta(i) = betax * Bloc(el, i) + betay * Cloc(el, i);
        end
        for i = 1:3
            for j = 1:3
                Sloc(i, j) = (delta(el) / (4 * Area(el))) * gradDotBeta(i) * gradDotBeta(j);
            end
        end

        % Assembling
        for i = 1:3
            iglob = triang(el, i);
            for j = 1:3
                jglob = triang(el, j);
                stabMat(iglob, jglob) = stabMat(iglob, jglob) + Sloc(i, j);
            end
        end
    end
end

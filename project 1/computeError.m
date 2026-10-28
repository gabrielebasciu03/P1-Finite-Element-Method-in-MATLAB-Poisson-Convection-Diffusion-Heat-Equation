function [err_h, Ai] = computeError(Nelem,Nodes,triang,coord,Area,uh,ureal)
% COMPUTEERROR  Estimate the L2 norm of the error via the trapezoidal
% (vertex/nodal quadrature) rule:
%
%   eps_h = ( sum_i (u_h(x_i) - u(x_i))^2 * A_i )^{1/2}
%
% where A_i is the nodal area: one third of the area of every triangle
% sharing node i (the same weights used for the load-vector assembly).

Ai = zeros(Nodes,1);
for el = 1:Nelem
    nodes = triang(el,:);
    Ai(nodes) = Ai(nodes) + Area(el)/3;
end

uex = ureal(coord(:,1), coord(:,2));
err_h = sqrt( sum( (uh - uex).^2 .* Ai ) );
end

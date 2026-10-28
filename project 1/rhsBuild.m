
% RHSBUILD  Assembles the rhs vector according to the project formulas:
% - forcing term f: bi = Ai * f(xi,yi)
% - Neumann contribution: bi = (1/2) * q(yi) * li

  function rhs = rhsBuild(Nodes, coord, Area, triang, force, NeuNod)
    rhs = zeros(Nodes, 1);

    % --- 1. % Contribution of the forcing term ---
    nodal_area = calculateNodalArea(Nodes, triang, Area);
    x = coord(:,1);
    y = coord(:,2);
    rhs = rhs + nodal_area .* force(x, y);

   

    % --- 2. Contribution Neumann ---
    li = calculateBoundaryLength(Nodes, triang, coord);

   for k = 1:length(NeuNod)
    idx = NeuNod(k);
    xi = coord(idx,1);
    yi = coord(idx,2);

    if abs(xi - 1) < 1e-8
        q = -pi * cos((pi/2)*yi);
    elseif abs(xi + 1) < 1e-8
        q =  pi * cos((pi/2)*yi);
    else
        q = 0;
    end

   

    
    rhs(idx) = rhs(idx) + q * li(idx);   
    end
end
%end
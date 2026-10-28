function rhs = rhsBuild2(Nodes, coord, Area, triang, force)
% RHSBUILD2  Builds the right-hand side vector for Project 2.
%     The "force" parameter is optional: if omitted, rhs = 0 everywhere
%     (the case required by Project 2). If instead you want to verify the
%     code against the Poisson solver of Project 1 (as
%     suggested by the text: "compare its solution with the solution
%     obtained with the available Poisson solver"), you can pass the
%     same force(x,y) used in Project 1.

    rhs = zeros(Nodes, 1);

    if nargin < 5 || isempty(force)
        return
    end

    nodal_area = calculateNodalArea(Nodes, triang, Area);
    x = coord(:, 1);
    y = coord(:, 2);
    rhs = rhs + nodal_area .* force(x, y);
end

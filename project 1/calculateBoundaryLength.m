function li = calculateBoundaryLength(Nodes, triang, coord)
% CALCULATEBOUNDARYLENGTH  Computes, for each node i, the quantity li,
% defined as the sum of half the length of the boundary edges
% (on Gamma_N: x = +-1) that share node i.

    li = zeros(Nodes, 1);
    tol = 1e-8;

% all triangle edges (each boundary edge appears only once,
% since it is shared by only one triangle)
    edges = [triang(:,[1 2]); triang(:,[2 3]); triang(:,[3 1])];

    for e = 1:size(edges,1)
        n1 = edges(e,1);
        n2 = edges(e,2);
        x1 = coord(n1,1);
        x2 = coord(n2,1);

        on_right = (abs(x1-1) < tol) && (abs(x2-1) < tol);   % x=1
        on_left  = (abs(x1+1) < tol) && (abs(x2+1) < tol);   % x=-1

        if on_right || on_left
            len = norm(coord(n1,:) - coord(n2,:));
            li(n1) = li(n1) + len/2;
            li(n2) = li(n2) + len/2;
        end
    end
end
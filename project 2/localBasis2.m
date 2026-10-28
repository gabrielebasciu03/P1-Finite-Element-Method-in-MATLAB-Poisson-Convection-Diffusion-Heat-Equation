function [Bloc, Cloc, Area] = localBasis2(Nelem, triang, coord)
% LOCALBASIS2  Computes the local geometric coefficients b_j, c_j of the
% P1 basis functions and the area of each triangle in the mesh.

    Bloc = zeros(Nelem, 3);
    Cloc = zeros(Nelem, 3);
    Area = zeros(Nelem, 1);

    for el = 1:Nelem
        nodes = triang(el, :);
        x = coord(nodes, 1);
        y = coord(nodes, 2);

        % Matrix for the determinant
        DetMat = [1, x(1), y(1); 1, x(2), y(2); 1, x(3), y(3)];
        Area(el) = 0.5 * abs(det(DetMat));

        % Coefficients bj = ym - yi, cj = xi - xm 
        Bloc(el, 1) = y(2) - y(3);
        Bloc(el, 2) = y(3) - y(1);
        Bloc(el, 3) = y(1) - y(2);

        Cloc(el, 1) = x(3) - x(2);
        Cloc(el, 2) = x(1) - x(3);
        Cloc(el, 3) = x(2) - x(1);
    end
end

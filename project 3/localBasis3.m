function [Bloc, Cloc, Area] = localBasis3(Nelem, triang, coord)
% LOCALBASIS3, P1 basis coefficients b,c and triangle areas (same as Project 1).
%   phi_k(x,y) = (a_k + b_k*x + c_k*y) / (2*Delta)

    Bloc = zeros(Nelem,3);
    Cloc = zeros(Nelem,3);
    Area = zeros(Nelem,1);

    for el = 1:Nelem
        nodes = triang(el,:);
        xy    = coord(nodes,:);

        Area(el) = 0.5 * abs( det([1 xy(1,:); 1 xy(2,:); 1 xy(3,:)]) );

        for k = 1:3
            k1 = mod(k,  3) + 1;
            k2 = mod(k+1,3) + 1;
            Bloc(el,k) = xy(k1,2) - xy(k2,2);
            Cloc(el,k) = xy(k2,1) - xy(k1,1);
        end
    end
end

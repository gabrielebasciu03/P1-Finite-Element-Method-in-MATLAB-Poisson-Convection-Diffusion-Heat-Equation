function Stiff = stiffBuild3(Nelem, Nodes, triang, Bloc, Cloc, Area)
% STIFFBUILD3, global stiffness matrix:
%   A_ij = (b_i*b_j + c_i*c_j) / (4*Delta)

    I_idx = zeros(Nelem*9,1);
    J_idx = zeros(Nelem*9,1);
    V     = zeros(Nelem*9,1);
    cnt   = 0;

    for el = 1:Nelem
        nodes = triang(el,:);
        Aloc  = zeros(3,3);
        for i = 1:3
            for j = 1:3
                Aloc(i,j) = (Bloc(el,i)*Bloc(el,j) + Cloc(el,i)*Cloc(el,j)) ...
                             / (4*Area(el));
            end
        end
        for i = 1:3
            for j = 1:3
                cnt = cnt + 1;
                I_idx(cnt) = nodes(i);
                J_idx(cnt) = nodes(j);
                V(cnt)     = Aloc(i,j);
            end
        end
    end

    Stiff = sparse(I_idx(1:cnt), J_idx(1:cnt), V(1:cnt), Nodes, Nodes);
end
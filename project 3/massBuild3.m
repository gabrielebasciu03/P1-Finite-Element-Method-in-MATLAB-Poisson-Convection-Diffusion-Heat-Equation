function Mass = massBuild3(Nelem, Nodes, triang, Area)
% MASSBUILD3, global consistent mass matrix for P1 elements:
%   Mloc = Delta/12 * [2 1 1; 1 2 1; 1 1 2]

    I_idx = zeros(Nelem*9,1);
    J_idx = zeros(Nelem*9,1);
    V     = zeros(Nelem*9,1);
    cnt   = 0;

    Mref = [2 1 1; 1 2 1; 1 1 2];

    for el = 1:Nelem
        nodes = triang(el,:);
        Mloc  = (Area(el)/12) * Mref;
        for i = 1:3
            for j = 1:3
                cnt = cnt + 1;
                I_idx(cnt) = nodes(i);
                J_idx(cnt) = nodes(j);
                V(cnt)     = Mloc(i,j);
            end
        end
    end

    Mass = sparse(I_idx(1:cnt), J_idx(1:cnt), V(1:cnt), Nodes, Nodes);
end
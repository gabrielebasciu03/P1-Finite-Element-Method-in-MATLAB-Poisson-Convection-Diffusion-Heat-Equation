function M = imposeBC3(Nodes, DirNod, M)
% IMPOSEBC3, impose homogeneous Dirichlet BC on a matrix by the penalty method
% row/col zeroed, diagonal = 1/eps.


    eps_pen = 1e-15;
    for k = 1:length(DirNod)
        nd = DirNod(k);
        M(nd,:) = 0;
        M(:,nd) = 0;
        M(nd,nd) = 1/eps_pen;
    end
    M = sparse(M);
end
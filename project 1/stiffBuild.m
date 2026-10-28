function stiffMat = stiffBuild(Nelem, Nodes, triang, Bloc, Cloc, Area)
    stiffMat = sparse(Nodes, Nodes);
    
    for el = 1:Nelem
        % Local stiffness matrix 3x3
        Sloc = zeros(3, 3);
        for i = 1:3
            for j = 1:3
                Sloc(i,j) = (1 / (4 * Area(el))) * (Bloc(el,i)*Bloc(el,j) + Cloc(el,i)*Cloc(el,j));
            end
            if Area(el) <= 0
               error('Not valid Area in %d. check the la mesh!', el);
          end 
    end
        
        % Assembing
        for i = 1:3
            iglob = triang(el, i);
            for j = 1:3
                jglob = triang(el, j);
                stiffMat(iglob, jglob) = stiffMat(iglob, jglob) + Sloc(i, j);
            end
        end
    end
end
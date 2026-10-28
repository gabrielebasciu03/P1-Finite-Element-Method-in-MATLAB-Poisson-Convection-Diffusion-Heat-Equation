function nodeArea = nodeAreaAssemble(Nelem, Nodes, triang, Area)
% NODEAREAASSEMBLE, lumped nodal areas, used only to compute the L2 spatial error norm: A_i = sum_{T ni} Area(T)/3

    nodeArea = zeros(Nodes,1);
    for el = 1:Nelem
        nodes = triang(el,:);
        for k = 1:3
            nodeArea(nodes(k)) = nodeArea(nodes(k)) + Area(el)/3;
        end
    end
end
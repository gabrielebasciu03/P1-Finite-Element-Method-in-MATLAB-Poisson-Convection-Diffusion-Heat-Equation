function nodal_area = calculateNodalArea(Nodes, triang, Area)
% CALCULATENODALAREA Auxiliary function to sum 1/3 of the area of the triangles incident on the node

    nodal_area = zeros(Nodes, 1);
    Nelem = size(triang, 1);
    for el = 1:Nelem
        for i = 1:3
            idx = triang(el, i);
            nodal_area(idx) = nodal_area(idx) + Area(el)/3;
        end
    end
end


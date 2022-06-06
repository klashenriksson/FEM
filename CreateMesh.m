function [p,e,t, bdry_nodes] = CreateMesh(g, h, polygrad, bdry_checker)
    [p,e,t] = initmesh(g, 'hmax', h);
    if polygrad == 2
        [p,t] = ChangeP1toP2Mesh(p,t);
    end

    bdry_nodes = [];
    for i = 1:size(p, 2)
        P = p(:,i);
        if bdry_checker(P(1), P(2)) == 1
            bdry_nodes(end+1) = i;
        end
    end

    t(end,:)=[];
end
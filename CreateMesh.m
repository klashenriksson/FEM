function [p,e,t] = CreateMesh(g, h, polygrad)
    [p,e,t] = initmesh(g, 'hmax', h);
    if polygrad == 2
        [p,t] = ChangeP1toP2Mesh(p,t);
    end

    t(end,:)=[];
end
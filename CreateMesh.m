function [p,e,t, bdry_nodes] = CreateMesh(g, h, polygrad, bdry_checker)
    [p,e,t] = initmesh(g, 'hmax', h);
    [p,e,t,bdry_nodes] = ConfigMesh(p,e,t,polygrad,bdry_checker);
end
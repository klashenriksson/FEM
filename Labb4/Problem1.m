%clear
%close all
rho = 1e-5;
gD = @(x,y) [0,0];
gN = @(x,y) [0,0];
f = @(x,y) rho*[0,-9.82]';
nu = 0.3;
ED = 1;
polygrad = 1;

g = Rectg(0,0,5,0.25); % unit square
h = 0.1;
[p,e,t] = CreateMesh(g, h, polygrad);
e = e(1:2,:);
edges = unique(e);
Dirichlet_boundary_nodes = [];
for E = 1:length(edges)
    if p(1,edges(E)) == 0
        Dirichlet_boundary_nodes(end+1) = edges(E);
    end
end


[U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p,e,t,f,gD,gN,ED,nu,polygrad,Dirichlet_boundary_nodes);

figure;
plotElasticity(p,t,U,D)
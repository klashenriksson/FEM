clear
close all
gD = @(x,y) [0;0];
gN = @(x,y) [0;0];
nu = 0.3;
ED = 1;
[mu, lambda] = Enu2Lame(ED, nu);
f = @(x,y) [(lambda + mu)*(1-2*x)*(1-2*y); -2*mu*y*(1-y) - 2*(lambda + 2 * mu)*x*(1-x)];
u_anal = @(x,y) [0;-x.*(1-x).*y.*(1-y)];
polygrad = 2;

g = Rectg(0,0,1,1); % unit square
h = 0.025;
[p,e,t] = CreateMesh(g, h, polygrad);

e = e(1:2,:);
edges = unique(e);
Dirichlet_boundary_nodes = [];
for E = 1:length(edges)
    Dirichlet_boundary_nodes(end+1) = edges(E);
end

[U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p,e,t,f,gD,gN,ED,nu,polygrad,Dirichlet_boundary_nodes);
%%
nt = size(t, 2);
u_anal_vec = zeros(length(U), 1);
for K = 1:nt
    dofs = t(:,K);
    loc2glb = zeros(6*polygrad, 1);
    loc2glb(1:2:end-1) = 2 * dofs-1; loc2glb(2:2:end) = 2 * dofs;
    xtri = p(1, dofs);
    ytri = p(2, dofs);

    u_a = zeros(length(xtri)*2, 1);
    for i = 1:length(xtri)
        u_a_i = u_anal(xtri(i), ytri(i));
        u_a(i*2 - 1) = u_a_i(1);
        u_a(i*2) = u_a_i(2);
    end
    u_anal_vec(loc2glb) = u_a;
end

%%
figure;
plotElasticity(p,t,U,D)
figure;
plotElasticity(p, t, u_anal_vec, D);

%%
figure;
fsurf(@(x,y) -x.*(1-x).*y.*(1-y), [0 1 0 1])
figure;

tsurf = [t; ones(1, size(t, 2))];
pdesurf(p,tsurf,U(2:2:end));
figure;
pdesurf(p, tsurf, u_anal_vec(2:2:end));
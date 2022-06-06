clear
close all
gD = @(x,y) [0;0];
gN = @(x,y) [0;0];
nu = 0.3;
ED = 1;
[mu, lambda] = Enu2Lame(ED, nu);
f = @(x,y) [(lambda + mu)*(1-2*x)*(1-2*y); -2*mu*y*(1-y) - 2*(lambda + 2 * mu)*x*(1-x)];
u_anal = @(x,y) [0;-x.*(1-x).*y.*(1-y)];
diri_bdry_checker = @(x,y) 1 * (x == 1 || x == 0 || y == 0 || y == 1);
polygrad = 2;

g = Rectg(0,0,1,1); % unit square

%%
hs = [0.5,0.25,0.125,0.0625];
energy_p1 = zeros(length(hs), 1);
energy_p2 = zeros(length(hs), 1);

for i = 1:length(hs)
    h = hs(i);

    [p,e,t, bdry_nodes] = CreateMesh(g, h, 1, diri_bdry_checker);
    [U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p,e,t,f,gD,gN,ED,nu,1,bdry_nodes);
    energy_p1(i) = U'*b_tilde;

    [p,e,t, bdry_nodes] = CreateMesh(g, h, 2, diri_bdry_checker);
    [U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p,e,t,f,gD,gN,ED,nu,2,bdry_nodes);
    energy_p2(i) = U'*b_tilde;
end
%%
[p,e,t, bdry_nodes] = CreateMesh(g, 0.1, polygrad, diri_bdry_checker);
[U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p,e,t,f,gD,gN,ED,nu,polygrad,bdry_nodes);

%%
figure;
plot(hs, energy_p1);
hold on
plot(hs, energy_p2);
hold on
fplot(@(x) (lambda + 3*mu)/90, [min(hs), max(hs)]);
legend("P1", "P2", "Analytic");
xlabel("Mesh size h");
ylabel("Energy Norm");
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
fsurf(@(x,y) 0, [0 1 0 1])
figure;

tsurf = [t; ones(1, size(t, 2))];
pdesurf(p,tsurf,U(1:2:end-1));
figure;
pdesurf(p, tsurf, u_anal_vec(1:2:end-1));
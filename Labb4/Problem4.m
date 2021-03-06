clear
close all

Domain = load("mesh_prob4.mat");
p = Domain.p;
e = Domain.e;
t = Domain.t;
rho = 7700;
gD = @(x,y) [0;0];
gN = @(x,y) [0;-rho*9.82*5000] * (y > max(p(2,:))-0.001);
nu = 0.28;
ED = 207*1e9;
[mu, lambda] = Enu2Lame(ED, nu);
f = @(x,y) [0; -rho*9.82];
u_anal = @(x,y) [0;-x.*(1-x).*y.*(1-y)];
diri_bdry_checker = @(x,y) 1 * (x == min(p(1,:)) || x == max(p(1,:)));
polygrad = 2;

%%
[p,e,t, bdry_nodes] = ConfigMesh(p,e,t, polygrad, diri_bdry_checker);
[U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p,e,t,f,gD,gN,ED,nu,polygrad,bdry_nodes);



%%
figure;
plotElasticity(p,t,U,D);
xlabel("x")
ylabel("y")


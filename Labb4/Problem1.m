clear
%close all
rho = 1e-5;
gD = @(x,y) [0;0];
gN = @(x,y) [0;0];
f = @(x,y) rho*[0;-9.82];
boundary_checker = @(x,y) 1 * (x == 0);
nu = 0.3;
ED = 1;
polygrad = 2;

g = Rectg(0,0,5,0.75); % unit square
h = 0.5;
[p,e,t, bdry_nodes] = CreateMesh(g, h, polygrad,boundary_checker);
% e = e(1:2,:);
% edges = unique(e);
% Dirichlet_boundary_nodes = [];
% for E = 1:length(edges)
%     if p(1,edges(E)) == 0
%         Dirichlet_boundary_nodes(end+1) = edges(E);
%     end
% end


[U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p,e,t,f,gD,gN,ED,nu,polygrad,bdry_nodes);
%%
figure;
plotElasticity(p,t,U,D)
hold on
xlabel("x");
ylabel("y");
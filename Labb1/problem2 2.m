clear

%% b)
% choose u(x,y) = x^2/2
g_D = @(x,y) x^2/2;
a = @(x,y) 1;
f = @(x,y) -1;

h = 0.1;
[U, p, e, t, ~] = U_FEM(h, f, g_D, a);

%% c)
u_analytic = @(x,y) x.^2/2;

pdesurf(p,t,U)
hold on

X = linspace(0,1, 1/h);
Y = linspace(0,1,1/h)';
Z = (X.^2/2 + Y.*0);
surf(X,Y,Z);

%% d)
hs = [0.2, 0.1, 0.05, 0.025];
errors_l2 = zeros(length(hs), 1);
errors_energy = zeros(length(hs), 1);

for i = 1:length(hs)
    [U, p, ~, t, A] = U_FEM(hs(i), f, g_D, a);

    errors_l2(i) = L2Error2D(p,t,u_analytic, U);

    % |||u|||^2 = int_omega (x^2) dx = 1/3
    u_anal_energy = 1/3;
    errors_energy(i) = sqrt(u_anal_energy - U'*A*U);
end
%%
eoc_errors_l2 = EOC(errors_l2, hs);
eoc_errors_energy = EOC(errors_energy, hs);

figure;
loglog(hs, errors_l2);
hold on
loglog(hs, 2.*hs); % we expect h^2 convergence
legend('FEM solution', 'Analytical')
figure;
loglog(hs, errors_energy);
hold on
loglog(hs, 1.*hs); % we expect h^1 convergence
legend('FEM solution', 'Analytical')

disp('EOC L2 Error: ' + mean(eoc_errors_l2))
disp('EOC Energy Error: ' + mean(eoc_errors_energy))

function [U, p, e, t, A] = U_FEM(h, f, g_D, a)
g = Rectg(0,0,1,1); % unit square
[p,e,t] = initmesh(g,'hmax',h); % create mesh
e = e(1:2,:);
boundary_nodes = unique(e);
boundary_values = zeros(length(boundary_nodes), 1);
for i = 1:length(boundary_values)
    boundary_values(i) = g_D(p(1, boundary_nodes(i)), p(2, boundary_nodes(i)));
end

A = StiffnessAssembler2D(p, t, a);
b = LoadAssembler2D(p,t,f);

[A_00, b_1, remainingdofs] = LockDofs(A, b, boundary_nodes, boundary_values);

epsi_0 = A_00\b_1;
U = zeros(length(b), 1);
U(remainingdofs) = epsi_0;
U(boundary_nodes) = boundary_values;
end
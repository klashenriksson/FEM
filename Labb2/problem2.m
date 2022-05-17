clear
close all

gD = @(x,y) x^2/2;
gN = @(x,y) pi*cos(pi) .* (x == 1) - pi*cos(0) .* (x==0);
a = @(x,y) 1;
f = @(x,y) pi*pi*sin(pi*x);
kappa = @(x,y) (0) .* (x > 0);

h = 0.05;
[U, p, e, t, A_tilde, b_tilde] = U_FEM_robin(h, f, a, kappa, gD, gN);

[V,D]=eigs(A_tilde,10,'sa');
V_0 = V(:,1);

figure;
pdesurf(p, t, V_0);
xlabel('x')
ylabel('y')

%%

a = LoadAssembler2D(p, t, @(x,y) 1);

MEGAMATRIX = [A_tilde a; a' 0];
LONG_LOAD = [b_tilde; 1];

U_MAGIC = MEGAMATRIX\LONG_LOAD;

figure;
pdesurf(p,t, U_MAGIC(1:end-1));
xlabel('x')
ylabel('y')

%%

function [U, p, e, t, A_tilde, b_tilde] = U_FEM_robin(h, f, a, kappa, gD, gN)
g = Rectg(0,0,1,1); % unit square
[p,e,t] = initmesh(g,'hmax',h); % create mesh

e = e(1:2,:);

A = StiffnessAssembler2D(p, t, a);
b = LoadAssembler2D(p,t,f);

[R, r] = RobinAssembler2D(p, t, kappa, gD, gN);

A_tilde = A + R;
b_tilde = b + r;

U = A_tilde\b_tilde;

end
clear
close all

gD = @(x,y) x^2/2;
gN = @(x,y) 0;
a = @(x,y) 1;
f = @(x,y) -1;
kappa = @(x,y) (10^6) .* (x > -1);

h = 0.1;
[U, p, e, t, ~] = U_FEM_robin(h, f, a, kappa, gD, gN);

figure;
pdesurf(p, t, U)

function [U, p, e, t, A] = U_FEM_robin(h, f, a, kappa, gD, gN)
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





clear
close all

gD = @(x,y) 0.5*x.^2;
gN = @(x,y) 0;
a = @(x,y) 1;
f = @(x,y) -1 + 2*pi*pi*sin(pi.*x).*sin(pi.*y);
kappa = @(x,y) (10^6) .* (x > 1);

h = 0.05;
[U1, p, e, t, ~] = U_FEM_robin(h, f, a, kappa, gD, gN);

kappa = @(x,y) (10^6) .* (x > 0.5);
[U2, p, e, t, ~] = U_FEM_robin(h, f, a, kappa, gD, gN);

figure;
pdesurf(p, t, U1)
xlabel('x')
ylabel('y')
hold on
%pdesurf(p, t, U2)

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





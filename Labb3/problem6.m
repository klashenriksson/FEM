clear
close all

epsilon = @(x,y) 11;
u_anal = @(x,y) 1/(1-exp(-2/epsilon(x,y)))*(1-exp((y-1)/epsilon(x,y)));
f = @(x,y) 0;
g_D = @(x,y) u_anal(x,y);
b = @(x,y) [0;1];

h = 0.1;
[U, p, e, t, A] = ConvectionDiffusionSolver2D(h, f, g_D, epsilon, b);
pdesurf(p,t,U)

hold on
fsurf(u_anal, [0,1])
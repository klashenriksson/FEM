clear
close all

gD = @(x,y) x^2/2;
gN = @(x,y) -sin(1) .* (x == 1);
a = @(x,y) 1;
f = @(x,y) 1;
kappa = @(x,y) (0) .* (x > 0);

g = Rectg(0,0,1,1); % unit square
[p,e,t] = initmesh(g,'hmax',h); % create mesh

a = LoadAssembler2D(p, t, f);

U = zeros(size(p), 1);

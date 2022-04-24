function [p, t, Pf] = L2Projector2D(f, h)
g = Rectg(0,0,1,1); % unit square
[p,e,t] = initmesh(g,'hmax',h); % create mesh
M = MassAssembler2D(p,t); % assemble mass matrix
b = LoadAssembler2D(p,t,f); % assemble load vector
Pf = M\b; % solve linear system
pdesurf(p,t,Pf) % plot projection
end
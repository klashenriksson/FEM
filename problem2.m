clear
g_D = @(x,y) x;
a = @(x,y) 1;
f = @(x,y) 0;

h = 0.01;
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
U = zeros(length(b));
U(remainingdofs) = epsi_0;
U(boundary_nodes) = boundary_values;

pdesurf(p,t,U)

function [An, bn, remainingdofs] = LockDofs(A, b, dofs, values)
    n_total = length(b);
    remainingdofs = (1:n_total);
    remainingdofs(dofs) = [];

    An = A(remainingdofs, remainingdofs);
    bn = b(remainingdofs) - A(remainingdofs, dofs)*values;
end
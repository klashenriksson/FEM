function [U, p, e, t, A] = ConvectionDiffusionSolver2D(h, f, g_D, epsilon, b)
g = Rectg(0,0,1,1); % unit square
[p,e,t] = initmesh(g,'hmax',h); % create mesh
e = e(1:2,:);
boundary_nodes = unique(e);
boundary_values = zeros(length(boundary_nodes), 1);
for i = 1:length(boundary_values)
    boundary_values(i) = g_D(p(1, boundary_nodes(i)), p(2, boundary_nodes(i)));
end

A = StiffnessAssembler2D(p, t, epsilon);
CA = ConvectionAssembler2D(p,t, b);
b = LoadAssembler2D(p,t,f);

A = A + CA;
[A_00, b_1, remainingdofs] = LockDofs(A, b, boundary_nodes, boundary_values);

epsi_0 = A_00\b_1;
U = zeros(length(b), 1);
U(remainingdofs) = epsi_0;
U(boundary_nodes) = boundary_values;
end
function [U, p, e, t, A_tilde, b_tilde,D] = ElastiscityFEMSolver(p, e, t, f, gD, gN,E,nu,polygrad,Dirichlet_boundary_nodes)
    [mu,lambda] = Enu2Lame(E,nu);

    D = [lambda + 2*mu,lambda, 0;lambda, lambda + 2*mu,0;0,0,mu];
    nbf = 3*polygrad;
    boundary_values = zeros(length(Dirichlet_boundary_nodes) * 2, 1);
    for i = 1:length(Dirichlet_boundary_nodes)
        bdry_val = gD(p(1, Dirichlet_boundary_nodes(i)), p(2, Dirichlet_boundary_nodes(i)));
        boundary_values(i*2 - 1) = bdry_val(1);
        boundary_values(i*2) = bdry_val(2);
    end
    
    A = ElastiscityAssembler2D(p,t,D,nbf);
    b = ElastiscityLoadVector2D(p,t,f,nbf);
    bn = ElastiscityNeumanLoadVector(p,e,t,gN,polygrad);

    A_tilde = A;
    b_tilde = b+bn;

    diri_dofs = zeros(length(Dirichlet_boundary_nodes) * 2, 1);
    diri_dofs(1:2:end-1) = 2 * Dirichlet_boundary_nodes(:) - 1;
    diri_dofs(2:2:end) = 2 * Dirichlet_boundary_nodes(:);

    [A_00, b_1, remainingdofs] = LockDofs(A_tilde, b_tilde, diri_dofs, boundary_values);

    epsi_0 = A_00\b_1;
    U = zeros(length(b), 1);
    U(remainingdofs) = epsi_0;
    U(diri_dofs) = boundary_values;
end
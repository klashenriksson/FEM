function [U, p, e, t, A_tilde, b_tilde, res, l2_error] = NonLinFEMSolver(p, e, t, f, a, kappa, gD, gN, u_k, do_newton, u_anal)
    e = e(1:2,:);

    boundary_nodes = unique(e);
    boundary_values = zeros(length(boundary_nodes), 1);
    for i = 1:length(boundary_values)
        boundary_values(i) = gD(p(1, boundary_nodes(i)), p(2, boundary_nodes(i)));
    end
    
    l2_error = [];
    res = [];

    u_k(boundary_nodes) = boundary_values;
    b = LoadAssembler2D(p,t,f);

    tol = 1e-8;
    max_iters = 100;
    iter = 0;
    while(iter == 0 || (iter < max_iters && res(iter) > tol))
        
        A = NonLinStiffnessAssembler2D(p, t, a, u_k);
        nl = NonLinLoadVectorComp(A,u_k);

        if do_newton
            A = A + NewtonMatrix(p,t,a,u_k);
        end

        A_tilde = A;
        b_tilde = b - nl;

        [A_00, b_1, remainingdofs] = LockDofs(A_tilde, b_tilde, boundary_nodes, boundary_values);

        epsi_0 = A_00\b_1;
        dU = zeros(length(b), 1);
        dU(remainingdofs) = epsi_0;

        u_k = u_k + dU;

        Anew = NonLinStiffnessAssembler2D(p,t,a,u_k);
        [Anew, b_1, remainingdofs] = LockDofs(Anew, b, boundary_nodes, boundary_values);

        res(iter+1) = max(abs(b_1 - Anew*u_k(remainingdofs)));
        
        if exist('u_anal', 'var')
            l2_error(iter+1) = L2Error2D(p,t,u_anal,u_k);
        end
        iter = iter + 1;
    end

    U = u_k;
end


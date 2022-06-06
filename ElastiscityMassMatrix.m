function M = ElastiscityMassMatrix(p, t, rho, nbf)
    np = size(p,2);
    nt = size(t,2);
    M = sparse(np*2,np*2); % allocate stiffness matrix

    for K = 1:nt
        M_K = zeros(2*nbf, 2*nbf);
        dofs = t(:,K);
        loc2glb = zeros(2*nbf, 1);
        loc2glb(1:2:end-1) = 2 * dofs-1; loc2glb(2:2:end) = 2 * dofs;
        xtri = p(1, t(1:3,K));
        ytri = p(2, t(1:3,K));

        J = [
                xtri(2) - xtri(1) ,xtri(3) - xtri(1);
                ytri(2) - ytri(1) ,ytri(3) - ytri(1) 
            ];
        area = 0.5 * det(J) ; % triangle area
        [rspts, qwgts] = Gausspoints(4);

        for iq = 1:length(qwgts) % loop over quadrature points
            r = rspts(iq,1); % local quadrature point
            s = rspts(iq,2);
            dA = area * qwgts(iq) ; % quadrature weight

            if nbf == 3
                [S, dSdr, dSds] = P1shapes(r,s);
            else
                [S, dSdr, dSds] = P2shapes(r,s);
            end

            dSdX = J'\[dSdr';dSds'];
            dSdx = dSdX(1,:)'; dSdy = dSdX(2,:)';

            [phi, ~] = BasisDisplacementAndStrain([S';dSdx';dSdy']);
            M_K = M_K + rho * (phi' * phi) .* dA;
        end

        M(loc2glb, loc2glb) = M(loc2glb, loc2glb) + M_K;
    end
end
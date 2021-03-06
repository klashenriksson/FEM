function b = ElastiscityLoadVector2D(p,t,f, nbf)
np = size(p,2);
nt = size(t,2);
b = zeros(np*2,1);

for K = 1:nt
    dofs = t(:,K);
    b_k = zeros(2*nbf, 1);
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
        x = J(1 ,:) * [r;s] + xtri(1); % global quadrature point
        y = J(2 ,:) * [r;s] + ytri(1);
        dA = area * qwgts(iq) ; % quadrature weight

        if nbf == 3
            [S, dSdr, dSds] = P1shapes(r,s);
        else
            [S, dSdr, dSds] = P2shapes(r,s);
        end
        
        dSdX = J'\[dSdr';dSds'];
        dSdx = dSdX(1,:)'; dSdy = dSdX(2,:)';

        [phi, ~] = BasisDisplacementAndStrain([S';dSdx';dSdy']);
        phif = phi' * f(x,y);
        b_k = b_k + phif .* dA;
    end
    b(loc2glb,1) = b(loc2glb,1) + b_k;
end
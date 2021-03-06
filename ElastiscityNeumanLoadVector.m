function [bn] = ElastiscityNeumanLoadVector(p,e,t,gN,polygrad)
    ne = size(e,2);
    np = size(p,2);
    [pts, qwgts] = Gausspoints1D(4); % line quadrature points pts in [0,1]
    nbf = 3*polygrad;
    bn = zeros(2*np,1);
    for E = 1:ne
        nodes = e(:,E);
        triangle = GetTriangle(nodes,t);
        loc2glb = zeros(2*nbf, 1);
        loc2glb(1:2:end-1) = 2*triangle-1;
        loc2glb(2:2:end) = 2*triangle;
        xtri = p(1,triangle(1:3));
        ytri = p(2,triangle(1:3));
        J = [
                xtri(2) - xtri(1) ,xtri(3) - xtri(1);
                ytri(2) - ytri(1) ,ytri(3) - ytri(1) 
            ];
        J_inv = inv(J);
        x1 = p(1,nodes(1)); x2 = p(1,nodes(2));
        y1 = p(2,nodes(1)); y2 = p(2,nodes(2));
        
        len = sqrt((x2-x1)^2 + (y2-y1)^2); % edge length
        
        be = zeros(2*nbf,1);
        for iq = 1: length(qwgts) % loop over quadrature points
  
            x = x1*pts(iq) + x2*(1-pts(iq)); % global quadrature point
            y = y1*pts(iq) + y2*(1-pts(iq));
            r = J_inv(1,:)*[x-x1; y-y1];
            s = J_inv(2,:)*[x-x1; y-y1];
            dS = len * qwgts(iq) ; % quadrature weight

            if nbf == 3
                [S, dSdr, dSds] = P1shapes(r,s);
            else
                [S, dSdr, dSds] = P2shapes(r,s);
            end
            [phi,~] = BasisDisplacementAndStrain([S';dSdr';dSds']);
        
            be = be + dS * phi' * gN(x,y) ; % add to integral value
    
        end % loop over quadrature points
        bn(loc2glb) = bn(loc2glb) + be;
    end
end
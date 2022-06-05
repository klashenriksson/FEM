function [bn] = ElastiscityNeumanLoadVector(p,e,gN,BDAS)
    ne = size(e,1);
    [pts, qwgts] = Gausspoints1D(4); % line quadrature points pts in [0,1]
    nbf = 2;
    bn = zeros(2*ne,1);
    for E = 1:ne
        loc2glb = zeros(2*nbf, 1);
        loc2glb(1:2:end-1) = 2*e(:,E)-1;
        loc2glb(2:2:end) = 2*e(:,E);
        nodes = e(:,E);
        x1 = p(1,nodes(1)); x2 = p(1,nodes(2));
        y1 = p(2,nodes(1)); y2 = p(2,nodes(2));
    
        len = sqrt((x2-x1)^2 + (y2-y1)^2); % edge length
        
        be = zeros(2*nbf,1);
        for iq = 1: length(qwgts) % loop over quadrature points
  
            x = x1*pts(iq) + x2*(1-pts(iq)); % global quadrature point
            y = y1*pts(iq) + y2*(1-pts(iq));
            [phi,~] = BDAS(pts(iq,1),pts(iq,2));
        
            dA = len*qwgts(iq); % quadrature weight
        
            be = be + dA .* gN(x,y) .* phi ; % add to integral value
    
        end % loop over quadrature points
        bn(loc2glb) = bn(loc2glb) + be;
    end
end
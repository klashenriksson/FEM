function R = RobinMassMatrix2D(p, e, kappa)
np = size(p,2); % number of nodes
ne = size(e,2); % number of boundary edges
R = sparse(np,np); % allocate boundary matrix

for E = 1:ne

    loc2glb = e(1:2,E); % boundary nodes
    xe = p(1,loc2glb); % node x-coordinates
    ye = p(2,loc2glb); % node y-

    len = sqrt((xe(1)-xe(2))^2+(ye(1)-ye(2))^2);

    phiL = @ (x,y) 1 - sqrt((x-xe(1))^2 + (y-ye(1))^2) / len;
    phiR = @ (x,y) 1 - sqrt((x-xe(2))^2 + (y-ye(2))^2) / len;
    
    RE = zeros(2, 2);
    
    RE(1,1) = NumIntegEdge(xe, ye, @(x,y) kappa(x,y)*phiL(x,y)*phiL(x,y));
    RE(2,2) = NumIntegEdge(xe, ye, @(x,y) kappa(x,y)*phiR(x,y)*phiR(x,y));
    RE(1,2) = NumIntegEdge(xe, ye, @(x,y) kappa(x,y)*phiL(x,y)*phiR(x,y));
    RE(2,1) = RE(1,2);
    
    R(loc2glb,loc2glb) = R(loc2glb,loc2glb) + RE;
end

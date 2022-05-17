function r = RobinLoadVector2D(p, e, kappa, gD, gN)
np = size(p,2);
ne = size(e,2);
r = zeros(np,1);
for E = 1:ne

    loc2glb = e(1:2,E);
    xe = p(1,loc2glb);
    ye = p(2,loc2glb);
    
    len = sqrt((xe(1)-xe(2))^2+(ye(1)-ye(2))^2);

    phiL = @ (x,y) 1 - sqrt((x-xe(1))^2 + (y-ye(1))^2) / len;
    phiR = @ (x,y) 1 - sqrt((x-xe(2))^2 + (y-ye(2))^2) / len;
    
    rE = zeros(2,1);

    rE(1) = NumIntegEdge(xe, ye, @ (x,y) (kappa(x,y)*gD(x,y) + gN(x,y))*phiL(x,y));
    rE(2) = NumIntegEdge(xe, ye, @ (x,y) (kappa(x,y)*gD(x,y) + gN(x,y))*phiR(x,y));

    r(loc2glb) = r(loc2glb) + rE;
end

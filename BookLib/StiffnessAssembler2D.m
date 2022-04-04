function A = StiffnessAssembler2D(p,t,a)
np = size(p,2);
nt = size(t,2);
A = sparse(np,np); % allocate stiffness matrix
for K = 1:nt
  loc2glb = t(1:3,K); % local-to-global map
  x = p(1,loc2glb); % node x-coordinates
  y = p(2,loc2glb); % node y-
  [area,b,c] = HatGradients(x,y);
  AK = NumIntegTrig(x,y,@(x,y) a(x,y).*(b.*b + c.*c));
  A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + AK; % add element stiffnesses to A
end
end

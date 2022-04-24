function A = StiffnessAssembler2D(p,t,a)
np = size(p,2);
nt = size(t,2);
A = sparse(np,np); % allocate stiffness matrix
for K = 1:nt
  loc2glb = t(1:3,K); % local-to-global map
  x = p(1,loc2glb); % node x-coordinates
  y = p(2,loc2glb); % node y-
  [area,b,c] = HatGradients(x,y);
  % grad phi_i = [b(1), c(1)]
  % grad phi_i+1 = [b(2), c(2)]
  phi_1 = [b(1), c(1)];
  phi_2 = [b(2), c(2)];
  phi_3 = [b(3), c(3)];

  phis = [phi_1;phi_2;phi_3];

  AK = zeros(size(phis, 1), size(phis, 1));
  for i=1:size(phis, 1)
      for j = 1:size(phis, 1)
          AK(i,j) = NumIntegTrig(x,y, @(x,y) a(x,y).*dot(phis(i,:), phis(j,:)));
      end
  end

  A(loc2glb,loc2glb) = A(loc2glb,loc2glb) + AK ; % add element stiffnesses to A
end

end

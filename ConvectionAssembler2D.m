function CA = ConvectionAssembler2D(p,t,b_vec)
np = size(p,2);
nt = size(t,2);
CA = sparse(np,np); % allocate stiffness matrix
for K = 1:nt
  loc2glb = t(1:3,K); % local-to-global map
  x = p(1,loc2glb); % node x-coordinates
  y = p(2,loc2glb); % node y-
  [~,a,b,c] = HatCoefficients(x', y', 0, 0);
  % grad phi_i = [b(1), c(1)]
  % grad phi_i+1 = [b(2), c(2)]
  grad_phi_1 = [b(1), c(1)];
  grad_phi_2 = [b(2), c(2)];
  grad_phi_3 = [b(3), c(3)];

  phi_1 = @(x,y) a(1) + b(1)*x + c(1)*y;
  phi_2 = @(x,y) a(2) + b(2)*x + c(2)*y;
  phi_3 = @(x,y) a(3) + b(3)*x + c(3)*y;

  grad_phis = [
      grad_phi_1;
      grad_phi_2;
      grad_phi_3
  ];
  phis = {phi_1, phi_2, phi_3};

  CAK = zeros(size(grad_phis, 1), size(grad_phis, 1));
  for i=1:size(grad_phis, 1)
      phi = phis{i};
      for j = 1:size(grad_phis, 1)
          CAK(i,j) = NumIntegTrig(x,y, @(x,y) dot(b_vec(x,y), grad_phis(j,:)).*phi(x,y));
      end
  end

  CA(loc2glb,loc2glb) = CA(loc2glb,loc2glb) + CAK ; % add element stiffnesses to A
end

end
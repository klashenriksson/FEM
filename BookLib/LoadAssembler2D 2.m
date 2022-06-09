function b = LoadAssembler2D(p,t,f)
np = size(p,2);
nt = size(t,2);
b = zeros(np,1);
for K = 1:nt
  loc2glb = t(1:3,K);
  x = p(1,loc2glb);
  y = p(2,loc2glb);

  bK_1 = NumIntegTrig(x, y, f, @(r,s) 1 - r - s);
  bK_2 = NumIntegTrig(x, y, f, @(r,s) r);
  bK_3 = NumIntegTrig(x, y, f, @(r,s) s);
  bK = [bK_1;bK_2;bK_3];
  b(loc2glb) = b(loc2glb) + bK;
end

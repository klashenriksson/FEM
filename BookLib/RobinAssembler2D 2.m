function [R,r] = RobinAssembler2D(p, t, kappa, gD, gN)

TR = triangulation(t(1:3 ,:)', p(1:2 ,:)');
e = freeBoundary(TR)'; % extracting free boundary

R = RobinMassMatrix2D(p, e, kappa);
r = RobinLoadVector2D(p, e, kappa, gD, gN);

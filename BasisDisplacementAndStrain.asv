function [Phi, B] = BasisDisplacementAndStrain(Smat)
S = Smat(1,:);
dSdr = Smat(2, :);
dSds = Smat(3,:);

nbf = length(S);
Phi = zeros(2, 2 * nbf); Phi(1, 1:2:end) = S; Phi(2, 2:2:end) = S;
B = zeros(3, 2*nbf); B(1,1:2:end) = dSdr; B(2, 2:2:end) = dSds;
B(3, 1:2:end) = dSds; B(3,2:2:end) = dSdr;
end
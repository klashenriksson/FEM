function [Phi,B] = GetPhiB(r,s,polygrad)
if polygrad ==1
    [S, dSdr, dSds] = P1shapes(r,s);
else
    [S, dSdr, dSds] = P2shapes(r,s);
end
smat = [S';dSdr';dSds'];
[Phi,B] = BasisDisplacementAndStrain(smat);
end


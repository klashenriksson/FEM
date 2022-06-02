function M = ShapeMat(r,s,p)
if p == 1
    [S, dSdr, dSds] = P1shapes(r,s);
else
    [S, dSdr, dSds] = P2shapes(r,s);
end

M = [S; dSdr; dSds];
end
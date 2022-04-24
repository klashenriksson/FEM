function value = NumIntegEdge(xe, ye, f)
%TR = triangulation(t(1:3 ,:)', p(1:2 ,:)');
%e = freeBoundary(TR)'; % extracting free boundary
%ned = size(e,2); % number of edges

[pts, qwgts] = Gausspoints1D (4); % line quadrature points pts in [0,1]
value = 0;

x1 = xe(1); x2 = xe(2);
y1 = ye(1); y2 = ye(2);

len = sqrt((x2-x1)^2 + (y2-y1)^2); % edge length

for iq = 1: length(qwgts) % loop over quadrature points
    x = x1*pts(iq) + x2*(1-pts(iq)); % global quadrature point
    y = y1*pts(iq) + y2*(1-pts(iq));

    dA = len*qwgts(iq); % quadrature weight

    value = value + dA * f(x,y); % add to integral value

end % loop over quadrature points

end


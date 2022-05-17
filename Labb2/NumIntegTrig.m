function value = NumIntegTrig(xtri, ytri, f, g_rs)
    J = [
            xtri(2) - xtri(1) ,xtri(3) - xtri(1);
            ytri(2) - ytri(1) ,ytri(3) - ytri(1) 
        ];
    area = 0.5 * det(J) ; % triangle area
    [rspts,qwgts] = Gausspoints(4) ; % ref triangle quadrature points
    value = 0;
    for iq = 1:length(qwgts) % loop over quadrature points
        r = rspts(iq,1); % local quadrature point
        s = rspts(iq,2);
        x = J(1 ,:) * [r;s] + xtri(1); % global quadrature point
        y = J(2 ,:) * [r;s] + ytri(1);
        dA = area * qwgts(iq) ; % quadrature weight

        if exist('g_rs', 'var')
            value = value + dA * f(x,y) * g_rs(r,s);
        else
            value = value + dA * f(x,y); % add to integral value
        end
    end
end
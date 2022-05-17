function err = L2Error2D(p,t,f,Pf)
    nel = size(t,2); % number of elements
    err2 = 0; % square of L2 - error
    for k = 1: nel % loop over elements
        idx = t(1:3,k);
        xpoints = p(1,idx);
        ypoints = p(2,idx);
        Pf_vals = Pf(idx);
        A = [xpoints(1), ypoints(1), 1;
            xpoints(2), ypoints(2), 1;
            xpoints(3), ypoints(3), 1
        ];
        Pf_coef = A\Pf_vals;

        f_err2 = @(x,y) (Pf_coef(1)*x + Pf_coef(2)*y + Pf_coef(3) - f(x,y)).^2;
        err2 = err2 + NumIntegTrig(xpoints, ypoints, f_err2);
    end % loop over elements
    err = sqrt (err2) ; % L2 - error
end
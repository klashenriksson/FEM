function [l2_err, h1_err] = Err2D(p,t,u_anal, grad_u_anal, uk)
    nel = size(t,2); % number of elements
    l2_err2 = 0; % square of L2 - error
    l2_grad_err2 = 0; % square of L2 error for gradients
    for k = 1: nel % loop over elements
        idx = t(1:3,k);
        xpoints = p(1,idx);
        ypoints = p(2,idx);
        uk_vals = uk(idx);
        A = [xpoints(1), ypoints(1), 1;
            xpoints(2), ypoints(2), 1;
            xpoints(3), ypoints(3), 1
        ];
        uk_coef = A\uk_vals;

        l2_err2_k = @(x,y) (uk_coef(1)*x + uk_coef(2)*y + uk_coef(3) - u_anal(x,y)).^2;
        l2_err2 = l2_err2 + NumIntegTrig(xpoints, ypoints, l2_err2_k);

        l2_grad_err2_k = @(x,y) dot([uk_coef(1);uk_coef(2)] - grad_u_anal(x,y), [uk_coef(1);uk_coef(2)] - grad_u_anal(x,y));
        l2_grad_err2 = l2_grad_err2 + NumIntegTrig(xpoints, ypoints, l2_grad_err2_k);

    end % loop over elements
    h1_err = sqrt(l2_err2 + l2_grad_err2);
    l2_err = sqrt(l2_err2);
end
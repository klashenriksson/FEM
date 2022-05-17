function N = NewtonMatrix(p,t,a_anom,uk)
    np = size(p,2);
    nt = size(t,2);
    N = sparse(np,np); % allocate Newton matrix
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
    
        tiny = 1e-8;
        ap_anom = @(x) (a_anom(x+tiny)-a_anom(x-tiny))/(2*tiny);
    
        phis = @(x,y,i) a(i) + b(i).*x + c(i).*y;
    
        grad_phis = [grad_phi_1;grad_phi_2;grad_phi_3];
    
        u_grad = [uk(loc2glb(1))*b(1) + uk(loc2glb(2))*b(2) + uk(loc2glb(3))*b(3); ...
                  uk(loc2glb(1))*c(1) + uk(loc2glb(2))*c(2) + uk(loc2glb(3))*c(3)];
    
        u_anom = @(x,y) uk(loc2glb(1))*(a(1) + b(1)*x + c(1)*y) ...
                       +uk(loc2glb(2))*(a(2) + b(2)*x + c(2)*y) ...
                       +uk(loc2glb(3))*(a(3) + b(3)*x + c(3)*y);
    
        NK = zeros(size(grad_phis, 1), size(grad_phis, 1));
        for i=1:size(grad_phis, 1)
            for j = 1:size(grad_phis, 1)
                NK(i,j) = NumIntegTrig(x,y, @(x,y) ap_anom(u_anom(x,y)).*phis(x,y,j).*dot(grad_phis(i,:), u_grad));
            end
        end
   
        N(loc2glb,loc2glb) = N(loc2glb,loc2glb) + NK ; % add element stiffnesses to N
    end
end
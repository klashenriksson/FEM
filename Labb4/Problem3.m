clear
close all

rho = 1;
g = Rectg(0,0,5,1);
h = 0.1;
polygrad = 1;
bdry_checker = @(x,y) 0;
E = 1;
nu = 0.3;
[mu,lambda] = Enu2Lame(E,nu);
D = [lambda + 2*mu,lambda, 0;lambda, lambda + 2*mu,0;0,0,mu];

[p,e,t,bdry_nodes] = CreateMesh(g, h, polygrad, bdry_checker);
A = ElastiscityAssembler2D(p, t, D, 3*polygrad);
% bignum = 10*max(diag(A)); 
% for i = 1:size(p,2)
%     x = p(1,i); y = p(2, i);
%     if x < 0.01 || x > 4.999
%         A(2*i-1, 2*i-1) = bignum; A (2*i, 2*i) = bignum;
%     end
% end
M = ElastiscityMassMatrix(p,t,rho,3*polygrad);

A = 0.5*(A + A');
M = 0.5*(M + M');

opts.isreal = 1; opts.issym=1;
[phi, omega2] = eigs(A, M, 10, 1e-6*(E/rho), opts);
[eigVals, I] = sort(diag(omega2));
eigModes = phi(:,I);

for i = 1:5
    figure;
    fac = 1;
    plotElasticity(p,t,fac*eigModes(:,i), D);
    title(strcat({'Eigenmode '},num2str(i)));
    xlabel(strcat({'Eigen frequency \omega='},...
    num2str(sqrt(eigVals(i))),' [rad/s]'));
end
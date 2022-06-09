% gives points and weights for integration of [0,1]
function [pts,qwgts] = Gausspoints1D(order)
if rem(order,2) == 0
   order = order+1;
end
n = (order+1)/2; % number of points
beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % 3-term recurrence coeffs
T = diag(beta,1) + diag(beta,-1);     % Jacobi matrix
[V,D] = eig(T);                       % Eigenvalue decomposition
x = diag(D); [x,i] = sort(x);         % Legendre points
w = 2*V(1,i).^2;                      % Quadrature weights
xi = (1+x)/2; % 0<xi<1
pts = xi;
qwgts = w/2; % sum w to 1

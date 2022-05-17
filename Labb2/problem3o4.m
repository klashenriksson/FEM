clear
close all

u_anal = @(x,y) sin(pi.*x).*sin(pi.*y);

gD = @(x,y) 0;
gN = @(x,y) 0;
a = @(x) 1 + x^2;
kappa = @(x,y) (10^6) .* (x > 0);

f = @(x,y)2*pi^2*sin(pi*x)*sin(pi*y)*(sin(pi*x)^2*sin(pi*y)^2+1) ...
                   - 2*pi^2*cos(pi*y)^2*sin(pi*x)^3*sin(pi*y) ...
                   - 2*pi^2*cos(pi*x)^2*sin(pi*x)*sin(pi*y)^3;

h = 0.05;
g = Rectg(0,0,1,1); % unit square
[p,e,t] = initmesh(g,'hmax',h); % create mesh

np = size(p, 2);
u_k = 100*rand(np,1);

%u_k = u_anal(p(1,:), p(2,:)) + 0.05;

[U_p, p, e, t, ~, ~, res_p, l2_error_p] = NonLinFEMSolver(p, e, t, f, a, kappa, gD, gN, u_k, 0, u_anal);
[U_n, p, e, t, ~, ~, res_n, l2_error_n] = NonLinFEMSolver(p, e, t, f, a, kappa, gD, gN, u_k, 1, u_anal);

%%

figure;
plot(log(l2_error_p))
xlabel('Iterations')
ylabel('log(l2\_error)')
figure;
plot(log(l2_error_n))
xlabel('Iterations')
ylabel('log(l2\_error)')

figure;
plot(log(res_p))
xlabel('Iterations')
ylabel('log(residual)')

figure;
plot(log(res_n))
xlabel('Iterations')
ylabel('log(residual)')

figure;
plot(log(l2_error_p))
hold on
plot(log(l2_error_n))
xlabel('Iterations')
ylabel('log(l2\_error)')
legend('Picard', 'Newton')

figure;
plot(log(res_p))
hold on
plot(log(res_n))
xlabel('Iterations')
ylabel('log(residual)')
legend('Picard', 'Newton')
%set(gca,'FontSize',13);



clear
close all

f = @(x,y) -2*pi*(cos(2*pi*x));
g = @(x,y) x + y;
A = @(x,y) [2 + sin(2*pi*x) 0; 0 2 + cos(2*pi*x)];
b = @(x,y) [0;0];
h = 0.1;
u_anal = @(x,y) x + y;
u_grad_anal = @(x,y) [1;1];

quadpts = [1,2,3,4];
l2_errs = zeros(numel(quadpts),1);
h1_errs = zeros(numel(quadpts),1);
for i = 1:numel(quadpts)
    quadpt = quadpts(i);
    [U, p, e, t, ~, L2Err, H1Err] = ConvectionDiffusionSolver2D(h,f,g,A, b, u_anal, u_grad_anal, quadpt);
    figure;
    pdesurf(p,t,U);
    hold on
    fsurf(u_anal, [0,1])
    l2_errs(i) = L2Err;
    h1_errs(i) = H1Err;
end

figure;
plot(quadpts, l2_errs);
xlabel("Gausspoints");
ylabel("L2 Error");

figure;
plot(quadpts, h1_errs);
xlabel("Gausspoints");
ylabel("H1 Error");


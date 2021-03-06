clear
close all

epsilon = @(x,y) 1;
u_anal = @(x,y) 1/(1-exp(-2/epsilon(x,y)))*(1-exp((y-1)/epsilon(x,y)));
u_grad_anal = @(x,y) [
    0;
    1/((1-exp(-2/epsilon(x,y)))*epsilon(x,y))*-exp((y-1)/epsilon(x,y));
];

f = @(x,y) 0;
g_D = @(x,y) u_anal(x,y);
b = @(x,y) [0;1];
hs = [0.8, 0.4, 0.2, 0.1, 0.05];

%%
l2_errs = zeros(length(hs), 1);
h1_errs = zeros(length(hs), 1);
for i = 1:length(hs)
    h = hs(i);
    [U, p, e, t, A, l2_err, h1_err] = ConvectionDiffusionSolver2D(h, f, g_D, epsilon, b, u_anal, u_grad_anal);
    l2_errs(i) = l2_err;
    h1_errs(i) = h1_err;
end

l2_eoc = EOC(l2_errs, hs);
h1_eoc = EOC(h1_errs, hs);

pdesurf(p,t,U)

hold on
fsurf(u_anal, [0,1])

figure;
plot(hs(1:end-1), l2_eoc)
xlabel("h");
ylabel("L2 EOC");
figure;
plot(hs(1:end-1), h1_eoc)
xlabel("h");
ylabel("H1 EOC");

%%
epsilons = [1, 10^(-1), 10^(-2), 10^(-3), 10^(-4)];
l2_errs = zeros(length(hs), numel(epsilons));
h1_errs = zeros(length(hs), numel(epsilons));
for i = 1:numel(epsilons)
    epsilon = @(x,y) epsilons(i);
    u_anal = @(x,y) 1/(1-exp(-2/epsilon(x,y)))*(1-exp((y-1)/epsilon(x,y)));
    u_grad_anal = @(x,y) [
        0;
        1/((1-exp(-2/epsilon(x,y)))*epsilon(x,y))*-exp((y-1)/epsilon(x,y));
    ];
    g_D = @(x,y) u_anal(x,y);
    for j = 1:numel(hs)
        h = hs(j);
        [U, p, e, t, A, l2_err, h1_err] = ConvectionDiffusionSolver2D(h, f, g_D, epsilon, b, u_anal, u_grad_anal);
        l2_errs(j,i) = l2_err;
        h1_errs(j,i) = h1_err;
    end
    figure;
    pdesurf(p,t,U);
    hold on
    fsurf(@(x,y) u_anal(x,y) - 0.3, [0,1]);
    xlabel("X");
    ylabel("Y");
    zlabel("U");
    title("\epsilon = " + string(epsilons(i)));
    legend("U_h", "U_{analytic}")
    set(gca, 'FontSize', 14)
end

%%

figure;
l2_eocs = zeros(numel(epsilons), numel(hs(1:end-1)));
h1_eocs = zeros(numel(epsilons), numel(hs(1:end-1)));
for i = 1:numel(epsilons)
    l2_errs_i = l2_errs(:,i);
    h1_errs_i = h1_errs(:,i);

    l2_eoc = EOC(l2_errs_i, hs);
    h1_eoc = EOC(h1_errs_i, hs);

    plot(hs(1:end-1), l2_eoc);
    hold on

    l2_eocs(i,:) = l2_eoc;
    h1_eocs(i,:) = h1_eoc;
end
xlabel("h");
ylabel("L2 EOC");
legend(string(epsilons))
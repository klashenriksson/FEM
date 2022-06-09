clear
close all


%% a
f = @ (x,y) sin(2*pi*x)*cos(2*pi*y);
[p, t, Pf] = L2Projector2D(f, 0.1);

%% c
hs = [0.2, 0.1, 0.05, 0.025];
l2_errors = zeros(length(hs), 1);
for i = 1:length(hs)
    h = hs(i);
    [p,t,Pf] = L2Projector2D(f, h);
    l2_errors(i) = L2Error2D(p,t,f,Pf);
end

EOCs = EOC(l2_errors, hs);
loglog(hs, l2_errors)
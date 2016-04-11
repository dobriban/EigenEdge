%% Partial solution to HW2, from STAT 325
cd('C:\Git\EigenEdge\Experiments\Stat 325\HW 2')
addpath('..\..\..\Code')
%% two-point mixture
gamma  = 1;
eps = 1/2;
w = [1-eps, eps];
t = [1, 30];
tic
[grid,density, ~, ~, mass_at_0] = compute_esd_ode(t, w, gamma);
toc
%%
plot(grid, density)


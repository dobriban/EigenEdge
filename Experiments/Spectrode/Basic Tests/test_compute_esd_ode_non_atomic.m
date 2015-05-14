%Test the ODE methods for solving MP equation
%this is a simple script with a few plots
%%
%% Setting 1: 
t=1/2;
r = [1 2];
w = 1/2;
w_int = 1/2;
%% Setting 2: 
t=2;
r = [0.5 1];
w = 1/2;
w_int = 1/2;
%%
gamma = 0.5;
epsi  = 1e-6;
[grid,density] =  compute_esd_ode_non_atomic(t,w,gamma, r,w_int,epsi);

figure,
plot(grid, density/max(density),'r','LineWidth',4)
%% empirical eigenvalues
p = 1e3;
n = floor(p/gamma);
pop_eigs = random_draw_disc_unif(t,w,r,w_int, p);
[X] = 1/sqrt(n)*randn(n,p)*diag(sqrt(pop_eigs));
D = svd(X,'econ').^2;
[heights,locations] = hist(D,3*floor(sqrt(p)));
width = locations(2)-locations(1);
heights = heights / max(heights);

% plot
figure
hold on
bar(locations,heights);
plot(grid, density/max(density),'r','LineWidth',4)
legend( 'Empirical Eigenvalues','Theoretical Prediction','Location','best');

%%
filename = sprintf( './non_atomic_spec.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
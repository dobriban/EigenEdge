%% examples for using the EigenEdge package, specifically, the Atomic function
%% Quick example
t = [1 5];
w = [1 1]/2;
gamma = 1/2;
[grid, density] =  compute_esd_ode(t, w, gamma);
figure, plot(grid, density/max(density),'r','LineWidth',4)
xlabel('Eigenvalue')
ylabel('Density');
%% Long example:  mixture of uniform + arithmetic model

gap = 0.01;
t_min = 1;
t_max = 10;
K = 10;
[t,w]  = arithmetic_model(K, t_min,t_max,gap);
w = 1/2*w;
t = 1+ t;

r = [0.5 1.5];
w_int = 1/2;

epsi = 10^(-6);
num_gp = 5*1e1;
gamma = 0.01;
[grid, density,~, ~, mass_at_0] =  compute_esd_ode_non_atomic(t, w, gamma, r,w_int,epsi);

%empirical eigenvalues
rng(0);
n = 5*1e4;
p = floor(gamma*n);
pop_eigs = random_draw_disc_unif(t,w,r,w_int, p);
[X] = 1/sqrt(n)*randn(n,p)*diag(sqrt(pop_eigs));
D = svd(X,'econ').^2;
%
[heights,locations] = hist(D,10*floor(sqrt(p)));
width = locations(2)-locations(1);
heights = heights / max(heights);

%plot
figure
hold on
bar(locations,heights);
plot(grid, density/max(density),'r','LineWidth',4)
set(gca,'fontsize',14)
xlabel('Eigenvalues');
h = legend( 'Empirical Eigenvalues','Theoretical Prediction','Location','best');
set(h,'FontSize',14);
%%
filename = sprintf( './Illustration_mixture_3.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
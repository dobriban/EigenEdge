%Two examples you can't do well otherwise
% mixture of point mass + uniform

%% Compute density
gap = 0.01;
t_min = 1;
t_max = 10;
K = 10;
[t,w]  = comb_model(K, t_min,t_max,gap);
w = 1/2*w;
t = 1+ t;

r = [0.5 1.5];
w_int = 1/2;

epsi = 10^(-6);
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
set(gca,'fontsize',20)
xlabel('Eigenvalues');
h = legend( 'Empirical Eigenvalues','Theoretical Prediction','Location','best');
set(h,'FontSize',20);
%% save plot
filename = sprintf( './Illustration_mixture.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% grayscale
col = colormap(gray(10));
figure
hold on
h = bar(locations,heights);
%set(h,'color',col(1,:));
h = plot(grid, density/max(density),'LineWidth',5);
set(h,'color',col(7,:));
set(gca,'fontsize',20)
xlabel('Eigenvalues');
h = legend( 'Empirical Eigenvalues','Theoretical Prediction','Location','best');
set(h,'FontSize',20);

%%
filename = sprintf( './bw_Illustration_mixture.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% compute functionals
num_gp = 50;
gamma_array = linspace(1e-2,0.4,num_gp)';

means = zeros(num_gp,1);
medians = zeros(num_gp,1);
modes = zeros(num_gp,1);

for i=1:length(gamma_array)
  tic
  gamma = gamma_array(i);
  [grid, density,~, ~, mass_at_0] =  compute_esd_ode_non_atomic(t, w, gamma, r,w_int,epsi);
  means(i) = esd_moment_grid(grid,density,mass_at_0,@(x)x);
  medians(i) = esd_quantile_grid(grid,density,mass_at_0,1/2);
  modes(i) = esd_mode_grid(grid,density);
  fprintf('%d-th iteration, time = %f \n', i, toc);
end

%% plot functionals as a function of gamma
a = {'-','--',':','-.'};
figure
hold on
h = plot(gamma_array,means,'linewidth',2); set(h,'LineStyle',a{1});
h = plot(gamma_array,medians,'linewidth',2); set(h,'LineStyle',a{2});
h = plot(gamma_array,modes,'linewidth',2); set(h,'LineStyle',a{3});
xlabel('gamma');
set(gca,'fontsize',20)
h = legend('mean','median','mode');
set(h,'FontSize',20);
xlim([min(gamma_array) max(gamma_array)]);
%%
filename = sprintf( './examples_mean_mode_arithmetic.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

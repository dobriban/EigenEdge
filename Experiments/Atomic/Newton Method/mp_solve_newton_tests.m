%An example where Newton method works and one where it doesn't
%% doesn't work (1): multi-cluster model, convergence to wrong solution
%leading to 'cutting off' a cluster
%ouch, this gives an example where it doesn't calculate the right density
t =[1; 2; 4; 6; 9];
w = ones(length(t),1)/length(t);
gamma = 1/10;
[grid, density] = compute_esd_newton(t,w,gamma);

figure, hold on
hist(t,20)
ind = density<Inf;
h = plot(grid, density/max(density(ind)),'r');
set(h,'Linewidth',3)
xlabel('Eigenvalue');
set(gca,'fontsize',14)
h = legend( 'population SD','density','Location','best');
set(h,'FontSize',14);
%
filename = sprintf( './newton_method_bad_case.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% A closely related example that works
t =[1; 2; 4; 6];
w = ones(length(t),1)/length(t);
gamma = 1/10;
[grid, density] = compute_esd_newton(t,w,gamma);

figure, hold on
hist(t,20)
h = plot(grid, density/max(density),'r');
set(h,'Linewidth',3)
xlabel('Eigenvalue');
set(gca,'fontsize',14)
h = legend('population SD','density','Location','best');
set(h,'FontSize',14);

%
filename = sprintf( './newton_method_good_case.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
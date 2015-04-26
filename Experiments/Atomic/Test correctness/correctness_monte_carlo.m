%Test Monte Carlo method for approximating density
%% MC: Increase n,p, Fit kernel estimator to fixed number of iid samples
gamma = 1/2;
rng(0)
p_grid = 10.^(1:3);
empirical_density = zeros(max(p_grid),length(p_grid));
num_monte = 1e3;
grid = linspace(1e-3,5,max(p_grid))';
gamma_plus = (1+sqrt(gamma))^2;
gamma_minus = (1-sqrt(gamma))^2;
MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;
theor_density = MP_density(grid);
for j=1:length(p_grid)
    p = p_grid(j);
    n = floor(p/gamma);
    t = ones(p,1);
    empirical_density_array = zeros(length(grid),num_monte);
    tic
    for i=1:num_monte
        X = 1/sqrt(n) * randn(n,p)*diag(t);
        lambda = svd(X,'econ');
        empirical_density_array(:,i) = ksdensity(lambda.^2,grid,'kernel','epanechnikov');
        t1 = toc;
        fprintf('%d/%d MC iterations with p=%d, took %f sec.\n', i,num_monte,p,t1);
    end
    t1 = toc;
    fprintf('%d MC iterations with p=%d, took %f sec.\n', num_monte,p,t1);
    empirical_density(:,j) = mean(empirical_density_array')';
    figure, plot(grid,empirical_density(:,j),'linewidth',3,'color',rand(1,3)), hold on
    plot(grid,theor_density,'linewidth',3,'color',rand(1,3)), hold on
end

err= empirical_density-theor_density*ones(1,length(p_grid));
MSE = zeros(length(p_grid),1);
for i=1:length(MSE)
    MSE(i) = norm(err(:,i))/sqrt(p_grid(i));
end
ind = (grid>0.01)&(grid<1.1*(1+sqrt(gamma))^2*max(t));

%%
a = {'--','-.','-','-.'};
rng(0)
figure, hold on
for i=1:3
    h  = plot(grid(ind),log10(abs(err(ind,i))));
    set(h,'LineStyle',a{i},'linewidth',3,'color',rand(1,3));
end
set(gca,'fontsize',14)
h = legend('p=10','p=100','p=1000');
xlabel('x_i')
ylabel('\Delta(x_i,\epsilon)');
%figtitle('MC: Increase n,p, Fixed number of iid samples (1e3)');
set(h,'FontSize',14);
%% plot
filename = sprintf( './Correctness_MC_increasing_np.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% MC (1) : Fixed n,p, Average over increasing iid samples
p = 200;
gamma = 1/2;
t = ones(p,1);
n = floor(p/gamma);
rng(0)

p_grid = 10.^(1:3);
empirical_density = zeros(p,length(p_grid));

for j=1:length(p_grid)
    num_monte = p_grid(j);
    grid = linspace(1e-3,5,p)';
    empirical_density_array = zeros(length(grid),num_monte);
    tic
    for i=1:num_monte
        X = 1/sqrt(n) * randn(n,p)*diag(t);
        lambda = svd(X,'econ');
        empirical_density_array(:,i) = ksdensity(lambda.^2,grid,'kernel','epanechnikov');
    end
    t1 = toc;
    fprintf('%f MC iterations took %f sec.\n', num_monte,t1);
    empirical_density(:,j) = mean(empirical_density_array')';
    %figure, plot(grid,empirical_density(:,j),'linewidth',3,'color',rand(1,3))
end

gamma_plus = (1+sqrt(gamma))^2;
gamma_minus = (1-sqrt(gamma))^2;
MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;
theor_density = MP_density(grid);
err= empirical_density-theor_density*ones(1,length(p_grid));
MSE = zeros(length(p_grid),1);
for i=1:length(MSE)
    MSE(i) = norm(err(:,i))/sqrt(p);
end
ind = (grid>0.01)&(grid<1.1*(1+sqrt(gamma))^2*max(t));

%
figure, hold on
plot(grid(ind),log10(abs(err(ind,1))))
plot(grid(ind),log10(abs(err(ind,2))),'*r')
plot(grid(ind),log10(abs(err(ind,3))),'+g')
legend('10 samples','100 samples','1000 samples');
xlabel('Eigenvalue')
ylabel('Log10 of error: theoretical vs numerical density');
figtitle('MC: Fixed n,p, Average over iid samples');
%% plot
filename = sprintf( './Correctness_MC_fixed_np.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
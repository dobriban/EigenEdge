%% Partial solution to P6 of HW4, from STAT 325
%%Approximate sampling variance of the CDF of Wishart eigenvalues
%%Comparison  to scalar iid formula q(1-q) 
cd('C:\Git\EigenEdge\Experiments\Stat 325\HW 4')
addpath('..\..\..\Code')
%% get a sample of eigenvalues
gamma  = 1/2;
n = 500;
p = floor(gamma*n);
n_m = 1e4;
eigs = zeros(n_m,min(n,p));
tic
for i=1:n_m
    X = randn(n,p)/sqrt(n);
    s = svd(X,'econ');
    eigs(i,:) = s.^2;
end
toc
%%
save('hw4_p6_raw_data.mat')
%% What is the empirical variability in F_n(q), for a grid of quantiles q?
K = 100; eps = 5*1e-3;
t=1; w = 1; qs = linspace(eps,1-eps,K);
grid= esd_quantile(t, w, gamma, qs);
var_cdf = zeros(K,1);

for k=1:K
mean_cdf = mean(eigs' < grid(k));
var_cdf(k) = var(mean_cdf);
end
var_cdf_scaled = p^2*var_cdf;

%% plot variance, and a few candidate models
figure, hold on
h1 = plot(qs, var_cdf_scaled,'LineWidth',2);
xlabel('Quantile q')
ylabel('Var[F_n(q)]')

v_theo_0 = @(q) q.*(1-q);
h2 = plot(qs,v_theo_0(qs),'LineWidth',2);
h3 = plot(qs,4*v_theo_0(qs),'LineWidth',2);
h4 = plot(qs,v_theo_0(qs)+1/2,'LineWidth',2);
legend([h1 h2 h3 h4], {'Emp', 'q(1-q)', '4q(1-q)', 'q(1-q)+1/2'}, 'location', 'best')
set(gca,'FontSize',14);
%%
filename = sprintf( './variance_of_MP_quantiles_straw_gamma=%.2f.png',gamma);
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);    

%% plot variance, and other candidate models
figure, hold on
h1 = plot(qs, var_cdf_scaled,'LineWidth',2);
xlabel('Quantile q')
ylabel('Var[F_n(q)]')

v_theo_0 = @(q) q.*(1-q);
h2 = plot(qs,0.9*v_theo_0(qs).^(1/10),'LineWidth',2);
legend([h1 h2], {'Emp', '0.9[q(1-q)]^{0.1}'}, 'location', 'best')
set(gca,'FontSize',14);
%%
filename = sprintf( './variance_of_MP_quantiles_improved_gamma=%.2f.png',gamma);
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);    


%% plot residual
figure, hold on
v_theo_0 = @(q) q.*(1-q);
v_cand = 0.9*v_theo_0(qs).^(1/10);
h1 = plot(qs, var_cdf_scaled-v_cand','LineWidth',2);
xlabel('Quantile q')
ylabel('Var[F_n(q)]')
legend(h1, {'residual'}, 'location', 'best')
set(gca,'FontSize',14);
%%
filename = sprintf( './variance_of_MP_quantiles_residual_gamma=%.2f.png',gamma);
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);    


%Test support identification in comb model

%% set parameters of comb model
gap = 0.01;
t_min = 1/2;
t_max = 10;
K = 6;
[t,w]  = comb_model(K, t_min,t_max,gap);

%% test Spectrode
epsi_array = 10.^(-(1:6))';
gamma_array = 2.^(linspace(-5,-2,4))';
K_err = zeros(length(epsi_array),length(gamma_array));
l_err = inf(length(epsi_array),length(gamma_array));
u_err = inf(length(epsi_array),length(gamma_array));
%grid_width = inf(length(epsi_array),length(gamma_array));

for i=1:length(epsi_array)
  epsi = epsi_array(i);
  for j=1:length(gamma_array)
    gamma = gamma_array(j);
    [grid, density, m, v, mass_at_0, K_hat, l_hat, u_hat, x, f_hat] =  compute_esd_ode(t, w, gamma,epsi,'brent');
    tic
    %get gold standard: K, support boundary
    epsi_iter  = 1e-7;
    [~, ~, m_i] =  compute_esd_fp(t,w,gamma,epsi_iter,grid);
    thresh = epsi_iter;
    [K,l,u,ind_l,ind_u] = support_identification(grid(grid>0),m_i(grid>0),thresh);
    fprintf('Test %d/%d; Subtest %d/%d took %2.2f seconds\n',i, length(epsi_array), j, length(gamma_array), toc);
    
    K_err(i,j) = abs(K-K_hat);
    [l_err(i,j),u_err(i,j)] = error_on_grid(grid,K_hat, l_hat, u_hat, K,ind_l,ind_u);
  end
end

%% Visualization
%% plot K_err as a function of epsi
a = {'-','--',':','-.'};
figure
hold on
h = plot(-log10(epsi_array),K_err(:,1),'linewidth',2); set(h,'LineStyle',a{1});
h = plot(-log10(epsi_array),K_err(:,2),'linewidth',2); set(h,'LineStyle',a{2});
h = plot(-log10(epsi_array),K_err(:,3),'linewidth',2); set(h,'LineStyle',a{3});
h = plot(-log10(epsi_array),K_err(:,4),'linewidth',2); set(h,'LineStyle',a{4});

xlim([min(-log10(epsi_array)), max(-log10(epsi_array))])
xlabel('Accuracy');
ylabel('\Delta_K');
set(gca,'fontsize',20)
h = legend('\gamma = 1/2^5','\gamma = 1/2^4','\gamma = 1/2^3','\gamma = 1/2^2');
set(h,'FontSize',20);
%%
filename = sprintf( './identify_number_of_intervals_support.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% plot l_err as a function of epsi
a = {'-','--',':','-.'};
figure
subplot(1,2,1), hold on
h = plot(-log10(epsi_array),log10(l_err(:,1)),'linewidth',2); set(h,'LineStyle',a{1});
h = plot(-log10(epsi_array),log10(l_err(:,2)),'linewidth',2); set(h,'LineStyle',a{2});
h = plot(-log10(epsi_array),log10(l_err(:,3)),'linewidth',2); set(h,'LineStyle',a{3});
h = plot(-log10(epsi_array),log10(l_err(:,4)),'linewidth',2); set(h,'LineStyle',a{4});
xlabel('Accuracy');
ylabel('log_{10}(\Delta_l)');
set(gca,'fontsize',20)
h = legend('\gamma = 1/2^5','\gamma = 1/2^4','\gamma = 1/2^3','\gamma = 1/2^2');
set(h,'FontSize',20);
xlim([min(-log10(epsi_array)), max(-log10(epsi_array))]);

subplot(1,2,2), hold on
h = plot(-log10(epsi_array),log10(u_err(:,1)),'linewidth',2); set(h,'LineStyle',a{1});
h = plot(-log10(epsi_array),log10(u_err(:,2)),'linewidth',2); set(h,'LineStyle',a{2});
h = plot(-log10(epsi_array),log10(u_err(:,3)),'linewidth',2); set(h,'LineStyle',a{3});
h = plot(-log10(epsi_array),log10(u_err(:,4)),'linewidth',2); set(h,'LineStyle',a{4});
xlabel('Accuracy');
ylabel('log_{10}(\Delta_u)');
set(gca,'fontsize',20)
h = legend('\gamma = 1/2^5','\gamma = 1/2^4','\gamma = 1/2^3','\gamma = 1/2^2');
set(h,'FontSize',20);
xlim([min(-log10(epsi_array)), max(-log10(epsi_array))]);

%%
filename = sprintf( './identify_endpoints_support.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

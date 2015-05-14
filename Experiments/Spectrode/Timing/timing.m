%Time the methods for solving MP equation
%% Compute theoretical and numerical solutions for several accuracies
% Null case where all eigenvalues are equal to 1
p = 200;
gamma = 1/2;
K = 5;
epsi_array = linspace(1,K,K);
t = 1;
w = 1;
rng(0)
% times = zeros(length(epsi_array),2);
errors = zeros(length(epsi_array),2);
gamma_plus = (1+sqrt(gamma))^2;
gamma_minus = (1-sqrt(gamma))^2;
MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;

for i=1:length(epsi_array);
    epsilon = 10.^(-epsi_array(i));
    %Spectrode
    tic
    [grid,density_ode] =  compute_esd_ode(t,w,gamma,epsilon);
    times(i,1) = toc;
    fprintf('Test %d/%d: For accuracy %2.2e, Spectrode method took %2.2f seconds\n',i,length(epsi_array), epsilon, times(i,1));
    
    %Iter
    multiplier_num_iter = 1;
    tic
    [density_fp] = compute_esd_fp(t,w,gamma,epsilon,grid,multiplier_num_iter);
    times(i,2) = toc;
    fprintf('Test %d/%d: For accuracy %2.2e, Iter method took %2.2f seconds\n',i,length(epsi_array), epsilon, times(i,2));
    
    %store errors
    theor_density = MP_density(grid);
    errors(i,1)= mean(density_ode-theor_density);
    errors(i,2)= mean(density_fp-theor_density);
end

%to save time when editing the figure associated with these results, I have
%saved the data to disk. These can be re-loaded to generate a new figure.
%%
%save('timing_experiment_results_identity.mat');
%%
load('timing_experiment_results_identity.mat');

%%
rng(1);
figure
a = {'--','-.','-','-.'};

%figtitle('All eigenvalues equal to 1');
subplot(1,2,1), hold on
c1 = rand(1,3);
h = plot(log10(times(:,1)),'linewidth',3,'color', c1);
set(h,'LineStyle',a{1});

c2 = rand(1,3);
h = plot(log10(times(:,2)),'linewidth',3,'color',c2);
set(h,'LineStyle',a{2});
set(gca,'fontsize',20)
xlabel('-log_{10}(\epsilon)')
ylabel('t(\epsilon)');
xlim([1 length(times(:,1))]);

%xlabel('Correct Sig. Digits Requested');
%ylabel('Log10 Running time');
h = legend('Spectrode','FPA','location','Best');
set(h,'FontSize',20);


subplot(1,2,2), hold on
h = plot(log10(errors(:,1)),'linewidth',3,'color',c1);
set(h,'LineStyle',a{1});
h = plot(log10(errors(:,2)),'linewidth',3,'color',c2);
set(h,'LineStyle',a{2});
set(gca,'fontsize',20)

xlabel('-log_{10}(\epsilon)')
ylabel('\Delta(\epsilon)');
xlim([1 length(times(:,1))]);

%xlabel('Correct Sig Digits  Requested');
%ylabel('Correct Sig Digits  Produced');
%h = legend('Spectrode','FPA','location','Best');
%set(h,'FontSize',20);

%% save plot
filename = sprintf( './Timing_tests_numerical_error_null.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% Eigenvalues are an equal mixture of 1 and 8
p = 200;
gamma = 1/2;
K = 6;
epsi_array = linspace(1,K,K);
t = [1,8];
w = [1,1]/2;

tau = 8;
rng(0)
times = zeros(length(epsi_array),2);
errors = zeros(length(epsi_array),2);

for i=1:length(epsi_array);
    epsilon = 10.^(-epsi_array(i));
    %Spectrode
    tic
    [grid,density_ode] = compute_esd_ode(t,w,gamma,epsilon);
    times(i,1) = toc;
    fprintf('Test %d/%d: For accuracy %2.2e, Spectrode method took %2.2f seconds\n',i,length(epsi_array), epsilon, times(i,1));
    
    %Iter
    n = floor(p/gamma);
    multiplier_num_iter = 1;
    tic
    [density_fp] = compute_esd_fp(t,w,gamma,epsilon,grid,multiplier_num_iter);
    times(i,2) = toc;
    fprintf('Test %d/%d: For accuracy %2.2e, Fixed Point method took %2.2f seconds\n',i,length(epsi_array), epsilon, times(i,2));
    
    %store errors
    [m_theor] = SparseStieltjes2(grid,gamma,1/2,tau);
    theor_density = 1/pi*imag(m_theor);
    errors(i,1)= mean(density_ode-theor_density);
    errors(i,2)= mean(density_fp-theor_density);
end


%%
%save('timing_experiment_results_two-point.mat');
%%
load('timing_experiment_results_two-point.mat');
%%
rng(1);
figure
a = {'--','-.','-','-.'};

%figtitle('Equal mixture of 1 and 2');
subplot(1,2,1), hold on
c1 = rand(1,3);
h = plot(log10(times(:,1)),'linewidth',3,'color', c1);
set(h,'LineStyle',a{1});

c2 = rand(1,3);
h = plot(log10(times(:,2)),'linewidth',3,'color',c2);
set(h,'LineStyle',a{2});
set(gca,'fontsize',20)
xlabel('-log_{10}(\epsilon)')
ylabel('t(\epsilon)');
xlim([1 length(times(:,1))]);
%xlabel('Correct Sig. Digits Requested');
%ylabel('Log10 Running time');
%h = legend('Spectrode','FPA','location','Best');
%set(h,'FontSize',20);


subplot(1,2,2), hold on
h = plot(log10(errors(:,1)),'linewidth',3,'color',c1);
set(h,'LineStyle',a{1});
h = plot(log10(errors(:,2)),'linewidth',3,'color',c2);
set(h,'LineStyle',a{2});
set(gca,'fontsize',20)

xlabel('-log_{10}(\epsilon)')
ylabel('\Delta(\epsilon)');
xlim([1 length(times(:,1))]);

%xlabel('Correct Sig Digits  Requested');
%ylabel('Correct Sig Digits  Produced');
h = legend('Spectrode','FPA','location','NorthEast');
set(h,'FontSize',20);


%% save plot
filename = sprintf( './Timing_tests_numerical_error_two-point_mixture.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);


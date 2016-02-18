%% Unit tests for optimal LSS
%% Compare numerically with OMH LSS
%for ep_L = 3; ep_U = 6; expt takes ~70 sec
%a run with ep=8 takes ~1000 sec
gamma = 1/2; %aspect ratio gamma = p/n
%Null
n = 1e3;
p = floor(n*gamma);
t = ones(p,1); %null distrib
s_null = 1; %the null spike

%Alternative
k=5; %number of different spikes
mini = 1.1; maxi = 0.99*(1+sqrt(gamma));
spikes_arr = linspace(mini,maxi,k)'; %the alternative spike

%Array of epsilons
ep_L = 3;
ep_U = 7;
ep_arr = 10.^(-(ep_L:1:ep_U));
%ep_arr = 10.^(-(3:1:5));

err  = zeros(length(spikes_arr),length(ep_arr));
times  = zeros(length(spikes_arr),length(ep_arr));
print_iter=1;
timer=tic;
for j=1:length(ep_arr)
    epsi = ep_arr(j);
    for i=1:length(spikes_arr)
        if print_iter==1
            str = sprintf('Epsi %d out of %d. Spike %d out of %d.\n',j,length(ep_arr),i,length(spikes_arr));
            fprintf(str);
            toc(timer)
        end
        s_alt = spikes_arr(i);
        %LSS_comput_method = 'collocation';
        LSS_comput_method = 'diag_regularization';
        experimental_mode = 0;
        tic
        [grid,LSS] = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method,n,experimental_mode,epsi);
        times(i,j) = toc;
        
        %OMH LSS
        z_0 = @(x) (1+x)*(gamma+x)/x;
        sp = z_0(s_alt-1);
        OMH_lss = -log(max(sp-grid,0));
        ind_in = (grid<((1+sqrt(gamma))^2))&(grid>((1-sqrt(gamma))^2));
        OMH_lss(~ind_in) = 0;
        OMH_lss = OMH_lss+max(abs(OMH_lss(OMH_lss<0)));
        OMH_lss(~ind_in) = 0;
        OMH_lss = OMH_lss/max(abs(OMH_lss));
        
        err(i,j) = mean(abs(LSS(ind_in)-OMH_lss(ind_in)));
    end
end

%% Aggregate MAD and plot
a = {'-','--','-.',':'};     savefigs =1;
err_agg = mean(err);
rng(2);
figure, hold on
h1 = plot(-log10(ep_arr),log10(err_agg),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
legend(h1,{'log Err'},'location','Best')
xlabel('-log_{10} Epsi')
ylabel('log_{10} Err');
set(gca,'fontsize',20)
%% save figures
if savefigs==1
    filename = sprintf( './Img/MAD OMH LSS gamma = %.2f ep_U = %d.png',gamma,ep_U);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end

%% Plot MAD for individual spikes
rng(2);
figure, hold on
plot(-log10(ep_arr),log10(err'),'linewidth',4);
legend('1.1','1.25','1.4','1.55','1.7','location','Best')
xlabel('-log_{10} Epsi')
ylabel('log_{10} Err');
set(gca,'fontsize',20)
%% save figures
if savefigs==1
    filename = sprintf( './Img/MAD OMH LSS indiv gamma = %.2f ep_U = %d.png',gamma,ep_U);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end

%% Aggregate times and plot
times_agg = mean(times);
rng(2);
figure, hold on
h1 = plot(-log10(ep_arr),log10(times_agg),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
legend(h1,{'log time'},'location','Best')
xlabel('-log_{10} Epsi')
ylabel('log_{10} Sec');
set(gca,'fontsize',20)
%% save figures
if savefigs==1
    filename = sprintf( './Img/time OMH LSS gamma = %.2f ep_U = %d.png',gamma,ep_U);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end

%% Plot times for individual spikes
rng(2);
figure, hold on
plot(-log10(ep_arr),log10(times'),'linewidth',4);
legend('1.1','1.25','1.4','1.55','1.7','location','Best')
xlabel('-log_{10} Epsi')
ylabel('log_{10} Sec');
set(gca,'fontsize',20)
%% save figures
if savefigs==1
    filename = sprintf( './Img/time OMH LSS indiv gamma = %.2f ep_U = %d.png',gamma,ep_U);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end
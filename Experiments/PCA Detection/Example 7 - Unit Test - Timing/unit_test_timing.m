%Unit tests for optimal LSS
%Compare running time and accuracy of diagonal regularization and
%collocation
%% Compare numerically with OMH LSS+running time
a = {'-','--','-.',':'};  plotfigs =0;
gamma = 1/2; %aspect ratio gamma = p/n
%Null
n = 1e3;
p = floor(n*gamma);
t = ones(p,1); %null distrib
s_null = 1; %the null spike

%Alternative
%k=5; %number of different spikes
%mini = 1.1; maxi = 0.99*(1+sqrt(gamma));
k=1;
spikes_arr = 0.9*(1+sqrt(gamma));

%Array of epsilons
ep_L = 3;
ep_U = 6;
ep_arr = 10.^(-(ep_L:1:ep_U));
%ep_arr = 10.^(-(3:1:5));

err_c  = zeros(length(spikes_arr),length(ep_arr));
err_d  = zeros(length(spikes_arr),length(ep_arr));
time_c  = zeros(length(spikes_arr),length(ep_arr));
time_d  = zeros(length(spikes_arr),length(ep_arr));
print_iter=1;
timer=tic;
for j=1:length(ep_arr)
    epsi = ep_arr(j);
    for i=1:length(spikes_arr)
        if print_iter==1
            str = sprintf('Epsi %d out of %d. Spike %d out of %d.\n',j,length(ep_arr),i,length(spikes_arr));
            fprintf(str);
            toc(timer);
        end
        s_alt = spikes_arr(i);
        
        LSS_comput_method = 'collocation';
        experimental_mode = 0;
        tic
        [grid_c,LSS_c] = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method,n,experimental_mode,epsi);
        time_c(i,j) = toc;
        
        LSS_comput_method = 'diag_regularization';
        experimental_mode = 0;
        tic
        [grid_d,LSS_d] = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method,n,experimental_mode,epsi);
        time_d(i,j) = toc;
        
        %OMH LSS
        z_0 = @(x) (1+x)*(gamma+x)/x;
        sp = z_0(s_alt-1);
        
        %err of colloc
        OMH_lss =-log(max(sp-grid_c,0));
        ind_in_c = (grid_c<((1+sqrt(gamma))^2))&(grid_c>((1-sqrt(gamma))^2));
        OMH_lss(~ind_in_c) = 0;
        OMH_lss = OMH_lss+max(abs(OMH_lss(OMH_lss<0)));
        OMH_lss(~ind_in_c) = 0;
        OMH_lss = OMH_lss/max(abs(OMH_lss(OMH_lss<Inf)));
        err_c(i,j) = mean(abs(LSS_c(ind_in_c)-OMH_lss(ind_in_c)));
        
        if plotfigs==1
        figure, hold on
        h1= plot(grid_c,LSS_c,'linewidth',4,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
        z_0 = @(x) (1+x)*(gamma+x)/x;
        sp = z_0(s_alt-1);
        OMH_lss = -log(max(sp-grid_c,0));
        ind_in_c = (grid_c<((1+sqrt(gamma))^2))&(grid_c>((1-sqrt(gamma))^2));
        OMH_lss(~ind_in_c) = 0;
        OMH_lss = OMH_lss+max(abs(OMH_lss(OMH_lss<0)));
        OMH_lss(~ind_in_c) = 0;
        OMH_lss = OMH_lss/max(abs(OMH_lss));
        h3 = plot(grid_c,OMH_lss,'linewidth',4,'color',rand(1,3));
        set(h3,'LineStyle',a{3});
        legend([h1 h3],{'LSS_c','OMH LSS'},'location','Best')
        end
        
        %err of d_reg
        OMH_lss =-log(max(sp-grid_d,0));
        ind_in_d = (grid_d<((1+sqrt(gamma))^2))&(grid_d>((1-sqrt(gamma))^2));
        OMH_lss(~ind_in_d) = 0;
        OMH_lss = OMH_lss+max(abs(OMH_lss(OMH_lss<0)));
        OMH_lss(~ind_in_d) = 0;
        OMH_lss = OMH_lss/max(abs(OMH_lss));
        err_d(i,j) = mean(abs(LSS_d(ind_in_d)-OMH_lss(ind_in_d)));
    end
end
str = sprintf('Done.\n');
fprintf(str);
toc(timer);
%% Aggregate and plot - Mean abs deviation
 savefigs =1;
if k>1
    err_agg_c = mean(err_c);
    err_agg_d = mean(err_d);
else
    err_agg_c = err_c;
    err_agg_d = err_d;
end
rng(2);
figure, hold on
h1 = plot(-log10(ep_arr),log10(err_agg_c),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(-log10(ep_arr),log10(err_agg_d),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
legend([h1 h2],{'Colloc','Diag Reg'},'location','Best')
xlabel('-log_{10} Epsi')
ylabel('log_{10} Err');
set(gca,'fontsize',20)
%% save figures
if savefigs==1
    filename = sprintf( './Img/MAD OMH LSS colloc diag gamma = %.2f ep_U = %d.png',gamma,ep_U);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end
%% Aggregate and plot - Running time
if k>1
    time_agg_c = mean(time_c);
    time_agg_d = mean(time_d);
else
    time_agg_c = time_c;
    time_agg_d = time_d;
end

rng(2);
figure, hold on
h1 = plot(-log10(ep_arr),log10(time_agg_c),'linewidth',4,'color',rand(1,3));
set(h1,'LineStyle',a{1});
h2 = plot(-log10(ep_arr),log10(time_agg_d),'linewidth',4,'color',rand(1,3));
set(h2,'LineStyle',a{2});
legend([h1 h2],{'Colloc','Diag Reg'},'location','Best')
xlabel('-log_{10} Epsi')
ylabel('log_{10} Time');
set(gca,'fontsize',20)
%% save figures
if savefigs==1
    filename = sprintf( './Img/time OMH LSS colloc diag gamma = %.2f ep_U = %d.png',gamma,ep_U);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end
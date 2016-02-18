%% Compute and plot optimal LSS
%% Example - White Noise
a = {'-','--','-.',':'};     savefigs =0;
gamma = 1/2; %aspect ratio gamma = p/n

%Null
n = 1e3;
p = floor(n*gamma);
t = ones(p,1); %null distrib
s_null = 1; %the null spike

%Alternative (one spike)
%big experiment
k=20; %number of trials
mini = 1.1; maxi = 3;

spikes_arr = linspace(mini,maxi,k)'; %the alternative spike
num_spikes=1;

print_iter=1;
tic
for i=1:length(spikes_arr)
    if print_iter==1
        str = sprintf('Spike %d out of %d.\n',i,length(spikes_arr));
        fprintf(str);
        toc
    end
    s_alt = spikes_arr(i);
    %LSS_comput_method = 'collocation';
    LSS_comput_method = 'diag_regularization';
    experimental_mode = 0;
    [grid,LSS, density] = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method,n,experimental_mode);
    
    rng(2);
    figure, hold on
    xlim([min(grid),max(grid)]);
    %plot the portion of the LSS outside the bulk with separate style
    if (s_alt>1+sqrt(gamma))||(s_alt<abs(1-sqrt(gamma)))
        h1 = plot(grid,LSS,'linewidth',4,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
    else
        %indices inside the bulk support
        u = (1+sqrt(gamma))^2;
        l = abs(1-sqrt(gamma))^2;
        ind_in = min(grid<u,grid>l);
        r = rand(1,3);
        h1 = plot(grid(ind_in==1),LSS(ind_in==1),'linewidth',4,'color',r);
        set(h1,'LineStyle',a{1});
        hold on
        h = plot(grid(grid<=l),LSS(grid<=l),'linewidth',4,'color',r);
        set(h,'LineStyle',a{4});
        hold on
        h = plot(grid(grid>=u),LSS(grid>=u),'linewidth',4,'color',r);
        set(h,'LineStyle',a{4});
    end
    
    h2 = plot(grid,density/max(abs(density)),'linewidth',4,'color',rand(1,3));
    set(h2,'LineStyle',a{2});
    set(gca,'fontsize',20)
    legend([h1 h2],{'LSS','density'},'location','Best')
    xlabel('Eigenvalue')
    
    %if relevant, plot OMH LSS
    if (num_spikes==1)&&(s_alt<1+sqrt(gamma))&&(s_alt>1)
        z_0 = @(x) (1+x)*(gamma+x)/x;
        sp = z_0(s_alt-1);
        OMH_lss = -log(max(sp-grid,0));
        ind_in = (grid<((1+sqrt(gamma))^2))&(grid>((1-sqrt(gamma))^2));
        OMH_lss(~ind_in) = 0;
        OMH_lss = OMH_lss+max(abs(OMH_lss(OMH_lss<0)));
        OMH_lss(~ind_in) = 0;
        OMH_lss = OMH_lss/max(abs(OMH_lss));
        h3 = plot(grid,OMH_lss,'linewidth',4,'color',rand(1,3));
        set(h3,'LineStyle',a{3});
        legend([h1 h3 h2],{'LSS','OMH LSS','density'},'location','Best')
    end
    title(num2str(s_alt));
    %% save figures
    if savefigs==1
        filename = sprintf( './Img/LSS white noise gamma = %.2f spike = %.2f.png',gamma,s_alt);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
end


%Compute the optimal LSS for testing a null hypothesis
%H_0: spectrum = (1-eps)*H+eps*G_0, against
%H_1: spectrum = (1-eps)*H+eps*G_1
%for general H,G_i - mixtures of point masses
%save fo files
%% Example 1 - gamma = 1/2
savefigs =1; %set savefigs =1; to save figures
a = {'-','-.','--',':'}; 
gamma = 1/2; %aspect ratio gamma = p/n

%Null
n = 1e2;
p = floor(n*gamma);
J =2;
mini = 1; maxi = 3;
base_eigs = linspace(mini,maxi,J)'; %null distrib
t  =repmat(base_eigs,floor(p/J),1);
s_null = base_eigs; %the null spike

%Alternative (one spike)
k=15; %number of trials
mini = 0.1; maxi = 5;
spikes_arr = linspace(mini,maxi,k)'; %the alternative spike

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
    [grid,LSS,density,asy_effect_size,K_hat,l_hat,u_hat,upper_pt,ind_in,any_outside] = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method);
    
    rng(2);
    figure, hold on
    if any_outside==1
        h1 = plot(grid,LSS,'linewidth',4,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
    else
        %plot each bulk interval separately
        r = rand(1,3);
        for j=1:K_hat
            ind = min(grid>=l_hat(j),grid<=u_hat(j));
            h1 = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
            set(h1,'LineStyle',a{1});
        end
        %now plot intervals outside of the bulk
        %leftmost interval
        ind = grid<l_hat(1);
        h = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
        set(h,'LineStyle',a{4});
        %rightmost interval
        ind = grid>u_hat(K_hat);
        h = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
        set(h,'LineStyle',a{4});
        %inner intervals
        for j=1:K_hat-1
            ind = min(grid>u_hat(j),grid<l_hat(j+1));
            h = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
            set(h,'LineStyle',a{4});
        end
    end
    h2 = plot(grid,density/max(abs(density)),'linewidth',4,'color',rand(1,3));
    set(h2,'LineStyle',a{2});
    set(gca,'fontsize',20)
    legend([h1 h2],{'LSS','density'},'location','Best')
    xlabel('Eigenvalue')
    
    title(num2str(s_alt));
    %% save figures
    if savefigs==1
        filename = sprintf( './Img/LSS gamma = %.2f numClus = %d spike = %.2f.png',gamma,J,s_alt);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
end

%% Example 2
gamma = 1/10; %aspect ratio gamma = p/n

%Null
n = 1e2;
p = floor(n*gamma);
J =2;
mini = 1; maxi = 3;
base_eigs = linspace(mini,maxi,J)'; %null distrib
t  =repmat(base_eigs,floor(p/J),1);
s_null = base_eigs; %the null spike

%Alternative (one spike)
k=15; %number of trials
mini = 0.1; maxi = 5;
spikes_arr = linspace(mini,maxi,k)'; %the alternative spike

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
     [grid,LSS,density,asy_effect_size,K_hat,l_hat,u_hat,upper_pt,ind_in,any_outside] = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method);
    
  rng(2);
    figure, hold on
    if any_outside==1
        h1 = plot(grid,LSS,'linewidth',4,'color',rand(1,3));
        set(h1,'LineStyle',a{1});
    else
        %plot each bulk interval separately
        r = rand(1,3);
        for j=1:K_hat
            ind = min(grid>=l_hat(j),grid<=u_hat(j));
            h1 = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
            set(h1,'LineStyle',a{1});
        end
        %now plot intervals outside of the bulk
        %leftmost interval
        ind = grid<l_hat(1);
        h = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
        set(h,'LineStyle',a{4});
        %rightmost interval
        ind = grid>u_hat(K_hat);
        h = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
        set(h,'LineStyle',a{4});
        %inner intervals
        for j=1:K_hat-1
            ind = min(grid>u_hat(j),grid<l_hat(j+1));
            h = plot(grid(ind),LSS(ind),'linewidth',4,'color',r);
            set(h,'LineStyle',a{4});
        end
    end
    h2 = plot(grid,density/max(abs(density)),'linewidth',4,'color',rand(1,3));
    set(h2,'LineStyle',a{2});
    set(gca,'fontsize',20)
    legend([h1 h2],{'LSS','density'},'location','Best')
    xlabel('Eigenvalue')
    
    title(num2str(s_alt));
    %% save figures
    if savefigs==1
        filename = sprintf( './Img/LSS gamma = %.2f numClus = %d spike = %.2f.png',gamma,J,s_alt);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
end

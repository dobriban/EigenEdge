%Finite sample test of optimal LSS
%MC simulation in a Toeplitz null and one-spike alternative
%generate and save:
%1. histogram of population eigenvalues
%2. pointwise plot of LSS
%3. scree plot and histogram of sample eigenvalues
%4. distrib of LSS under null and alternative, for an MC simulation
%% Several plots: Toeplitz null and one-spike alternative
%% scree plot, LSS, and distrib of LSS under null and alternative.
a = {'-','--','-.',':'};       plotfigs=1;  savefigs =1;
gamma = 1/2; %aspect ratio gamma = p/n

%Null - Toeplitz
n=500;
p = floor(n*gamma);
rho = 0.5; %0.5. 0.7 or 0.85;
r = rho.^(0:1:p-1);
Sigma = toeplitz(r);
t = eig(Sigma); %null distrib
%for toeplitz: min  = (1-rho)/(1+rho); max = (1+rho)/(1-rho).
s_null = [1]; %the null spike

%Alternative (one spike)
% k=15; %number of trials
% mini = 0.1; maxi = 5;
k=1; %number of trials
%s1 = (1+rho)/(1-rho)*(1+sqrt(gamma))*0.3;
s1 = 3.5; %3.5 for rho=0.5; 4.1 for 0.6; 6.5 for 0.7 
spikes_arr =[s1]; %spikes

print_iter=1;
tic
for i=1:length(spikes_arr(:,1))
    if print_iter==1
        str = sprintf('Spike %d out of %d.\n',i,length(spikes_arr(:,1)));
        fprintf(str);
        toc
    end
    s_alt = spikes_arr(i,:)';
    %plot population eigenvalues
    if (plotfigs==1)
        [heights_n,locations_n] = hist([t;s_null],10*floor(sqrt(p)));
        heights_n = heights_n / max(heights_n);
        [heights_a,locations_a] = hist([t;s_alt],10*floor(sqrt(p)));
        heights_a = heights_a / max(heights_a);
        
        figure, hold on
        subplot(1,2,1)
        bar(locations_n,heights_n);
        set(gca,'fontsize',20)
        xlabel('Eigenvalues');
        title('Pop Null')
        subplot(1,2,2)
        bar(locations_a,heights_a);
        set(gca,'fontsize',20)
        xlabel('Eigenvalues');
        title('Pop Alt')
        
        if savefigs==1
            filename = sprintf( './Img/pop hist toeplitz gamma = %.2f rho = %.2f spike = %.2f h=%d.png',...
                gamma,rho,s_alt(1),length(s_alt));
            saveas(gcf, filename,'png');
            fprintf(['Saved Results to ' filename '\n']);
            close(gcf)
        end
    end
    %choose how to compute LSS
    LSS_comput_method = 'diag_regularization'; %fast
    %LSS_comput_method =  'collocation'; %slower, potentially more accurate
    [grid,LSS,density,asy_effect_size,K_hat,l_hat,u_hat,upper_pt,ind_in,any_outside]  = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method);
    power = normcdf(norminv(0.05)+asy_effect_size);
    
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
    xlabel('Eigenvalue')
set(gca,'fontsize',20)
    %% save pointwise plot of LSS
    if savefigs==1
        filename = sprintf( './Img/rho %.1f/LSS toeplitz gamma = %.2f rho = %.2f spike = %.2f h=%d.png',rho,gamma,rho,s_alt(1),length(s_alt));
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
    
    %% finite sample simulation
    nMonte = 200;
    tn =  zeros(nMonte,1); %store null and alt LSS test stats
    ta =  zeros(nMonte,1);
    trace_n =  zeros(nMonte,1); %store null and alt test stats - mock: trace
    trace_a =  zeros(nMonte,1);
    
    for l=1:nMonte
        if print_iter==1
            str = sprintf('nMonte %d out of %d.\n',l,nMonte);
            fprintf(str);
            toc
        end
        %null and alt simu: assumes s_null and s_alt have equal length
        X_0=randn(n,p+length(s_null));
        X_null = n^(-1/2)*X_0*diag([t; s_null].^(1/2));
        eig_null = svd(X_null).^2;
        [heights_n,locations_n] = hist(eig_null,10*floor(sqrt(p)));
        heights_n = heights_n / max(heights_n);
        
        %alt
        X_alt = n^(-1/2)*X_0*diag([t; s_alt].^(1/2));
        eig_alt = svd(X_alt).^2;
        [heights_a,locations_a] = hist(eig_alt,10*floor(sqrt(p)));
        heights_a = heights_a / max(heights_a);
        
        %plot histogram of eigenvalues for first MC iter
        if (plotfigs==1)&&(l==nMonte)
            figure, hold on
            subplot(1,2,1)
            bar(locations_n,heights_n);
            set(gca,'fontsize',20)
            xlabel('Eigenvalues');
            title('Emp Null')
            subplot(1,2,2)
            bar(locations_a,heights_a);
            set(gca,'fontsize',20)
            xlabel('Eigenvalues');
            title('Emp Alt')
            
            if savefigs==1
                filename = sprintf( './Img/rho %.1f/empir hist toeplitz gamma = %.2f rho = %.2f spike = %.2f h=%d nMonte=%d.png',...
                    rho,gamma,rho,s_alt(1),length(s_alt),nMonte);
                saveas(gcf, filename,'png');
                fprintf(['Saved Results to ' filename '\n']);
                close(gcf)
            end
            
            %scree plot  of eigenvalues
            figure, hold on
            subplot(1,2,1)
            nu = 10;
            ei = 100*eig_null(1:nu)/sum(eig_null);
            bar(1:nu,ei)
            set(gca,'fontsize',20)
            xlim([0 nu+1/2]);
            xlabel('PC')
            ylabel('% Var')
            ax = gca;
            ax.XTick = [0 5 10];
            subplot(1,2,2)
            ei = 100*eig_alt(1:nu)/sum(eig_alt);
            bar(1:nu,ei)
            set(gca,'fontsize',20)
            xlim([0 nu+1/2]);
            xlabel('PC')
            ylabel('% Var')
            ax = gca;
            ax.XTick = [0 5 10];
            
            if savefigs==1
                filename = sprintf( './Img/rho %.1f/scree toeplitz gamma = %.2f rho = %.2f spike = %.2f h=%d nMonte=%d.png',...
                    rho,gamma,rho,s_alt(1),length(s_alt),nMonte);
                saveas(gcf, filename,'png');
                fprintf(['Saved Results to ' filename '\n']);
                close(gcf)
            end
        end
        
        %now compute the value of the LSS
        LSS_null = interp1(grid,LSS,eig_null,'pchip');
        LSS_alt = interp1(grid,LSS,eig_alt,'pchip');
        
        tn(l) = sum(LSS_null);
        ta(l) = sum(LSS_alt);
        trace_n(l) =  sum(eig_null); %store null and alt test stats - mock: trace
        trace_a(l) =  sum(eig_alt);
    end
    %Hist of LSS values (aggregated over all MCs)
    figure, hold on
    subplot(1,2,1)
    hist((tn-mean(tn))/std(tn));
    set(gca,'fontsize',14)
    xlabel('Null');
    set(gca,'FontSize',20);
    subplot(1,2,2)
    ta_std=(ta-mean(tn))/std(tn);
    hist(ta_std);
    set(gca,'fontsize',14)
    xlabel('Alternative');
    set(gca,'FontSize',20);
    str = sprintf('m = %.2f s = %.2f',mean(ta_std),std(ta_std));
    title(str);
    if savefigs==1
        filename = sprintf( './Img/rho %.1f/LSS standardized hist toeplitz gamma = %.2f rho = %.2f spike = %.2f h=%d nMonte=%d.png',...
            rho,gamma,rho,s_alt(1),length(s_alt),nMonte);
        saveas(gcf, filename,'png');
        fprintf(['Saved Results to ' filename '\n']);
        close(gcf)
    end
    mean_LSS_alt = mean(ta_std);
    std_LSS_alt = std(ta_std);
    
    %Hist of trace values (for comparison)
    figure, hold on
    subplot(1,2,1)
    hist((trace_n-mean(trace_n))/std(trace_n));
    set(gca,'fontsize',14)
    xlabel('Null');
    set(gca,'FontSize',20);
    subplot(1,2,2)
    trace_a_std=(trace_a-mean(trace_n))/std(trace_n);
    hist(trace_a_std);
    set(gca,'fontsize',14)
    xlabel('Alternative');
    set(gca,'FontSize',20);
    str = sprintf('m = %.2f s = %.2f',mean(trace_a_std),std(trace_a_std));
    title(str);
    mean_trace_alt = mean(trace_n);
    std_trace_alt = std(trace_n);
    
end

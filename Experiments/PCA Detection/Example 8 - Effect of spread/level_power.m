%Finite sample test of optimal LSS:
%MC Simulation to ascertain power in a spiked Toeplitz model
%Compare with power of top eigenvalue based test

%Scale the white MP distribution
%% Toeplitz null and one-spike alternative
a = {'-','--','-.','-','-.'};     savefigs =1;     plotfigs=0;
gamma = 1/2; %aspect ratio gamma = p/n
alpha = 0.05;
%Null - Toeplitz
n = 5*1e2;%1e3;
p = floor(n*gamma);

%eig_model = 'ones';
eig_model = 'toeplitz';
%eig_model = 'scaled_ones';
switch eig_model
    case 'toeplitz'
        rho = 0.7;
        r = rho.^(0:1:p-1);
        Sigma = toeplitz(r);
        %toeplitz t: min  = (1-rho)/(1+rho); max = (1+rho)/(1-rho)
        t = eig(Sigma); %null distrib
        %Alternative (one spike)
        s1 = (1+rho)/(1-rho)*(1+sqrt(gamma));
        %note upper_pt = 6.3377
        s_null = 1; %the null spike
    case 'ones'
        t = ones(p,1); %null distrib
        rho = 0;
        s1 = 3;
        s_null = 1; %the null spike
    case 'scaled_ones'
        upper_pt = 6.3377; %the PT for Toeplitz, rho = 0.7
        sigma = upper_pt/(1+sqrt(gamma));
        t = sigma*ones(p,1); %null distrib
        rho = 0;
        
        %rho_pt: use the upper boundary as for Toeplitz
        rho_pt = 0.7;
        s1 = (1+rho_pt)/(1-rho_pt)*(1+sqrt(gamma));
        s_null = sigma; %the null spike
end


num_spikes=10; %5 or 20
spikes_arr =linspace(1.01*s_null,s1,num_spikes)'; %spikes
%spikes_arr = (1+sqrt(gamma))*1.1;
power_LSS= zeros(length(spikes_arr(:,1)),1);
power_eigmax = zeros(length(spikes_arr(:,1)),1);

%First get the null distribution/critical value for the top eigenvalue
nMonte = 5*1e2;
eig_max =  zeros(nMonte,1);

print_iter=1;
tic
for l=1:nMonte
    %null and alt simu: assumes s_null and s_alt have equal length
    X_0=randn(n,p+length(s_null));
    X_null = n^(-1/2)*X_0*diag([t; s_null].^(1/2));
    eig_null = svd(X_null).^2;
    eig_max(l) = eig_null(1);
end

%set empirical critical values for test statistics under the null
crit_val_eig_max= quantile(eig_max,1-alpha);

for i=1:length(spikes_arr(:,1))
    if print_iter==1
        str = sprintf('Spike %d out of %d.\n',i,length(spikes_arr(:,1)));
        fprintf(str);
        toc
    end
    s_alt = spikes_arr(i,:)';
    %choose how to compute LSS
    LSS_comput_method = 'diag_regularization'; %fast
    %LSS_comput_method =  'collocation'; %slower, potentially more accurate
    [grid,LSS, density,asy_effect_size,~,~,~,upper_pt] = optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method,p);
    asy_power = normcdf(norminv(0.05)+asy_effect_size);
    
    %% finite sample simulation
    tn =  zeros(nMonte,1); %store null and alt LSS test stats
    ta =  zeros(nMonte,1);
    eig_max =  zeros(nMonte,1);
    
    for l=1:nMonte
        if (print_iter==1)&&(mod(l,100)==1)
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
        eig_max(l) = eig_alt(1);
        [heights_a,locations_a] = hist(eig_alt,10*floor(sqrt(p)));
        heights_a = heights_a / max(heights_a);
        
        %now compute the value of the LSS
        LSS_null = interp1(grid,LSS,eig_null,'pchip');
        LSS_alt = interp1(grid,LSS,eig_alt,'pchip');
        tn(l) = sum(LSS_null);
        ta(l) = sum(LSS_alt);
    end
    %set empirical critical values for LSS test statistics under the null
    crit_val_LSS = quantile(tn,1-alpha);
    %need to deal with potentially equal values of the LSS
    truly_larger_n = mean(tn>crit_val_LSS); %a fraction at most alpha that is truly larger than crit_val_LSS
    equal_n = mean(tn==crit_val_LSS); %a fraction that is equal to crit_val
    %probability with which a randomized rule rejects when seeing crit_val
    prop=0;
    if equal_n>0
        prop = (alpha-truly_larger_n)/equal_n;
    end;
    %compute empirical power for LSS
    truly_larger_a = mean(ta>crit_val_LSS);
    equal_a  = mean(ta==crit_val_LSS);
    power_LSS(i) = truly_larger_a+equal_a*prop;
    
    %power of eig max
    power_eigmax(i) = mean(eig_max>crit_val_eig_max);
end

% Plot power
if length(spikes_arr(:,1))>1
    rng(2);
    figure, hold on
    h = plot(spikes_arr,power_eigmax,'linewidth',4,'color',rand(1,3));
    set(h,'LineStyle',a{1});
    h = plot(spikes_arr,power_LSS,'linewidth',4,'color',rand(1,3));
    set(h,'LineStyle',a{2});
    set(gca,'fontsize',20)
    h= legend('Eig_{max}','LSS','location','Best');
    xlabel('Spike')
    xlim([min(spikes_arr),max(spikes_arr)])
    
    %plot upper edge
    SP=upper_pt; %your point goes here
    x=[SP,SP];
    y=get(gca,'Ylim');
    plot(x,y,'linewidth',2)
end
% save figures
if savefigs==1
    filename = sprintf( './Img/power toeplitz gamma = %.2f rho = %.2f n = %d nMonte = %d.png', gamma,rho,n,nMonte);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
end


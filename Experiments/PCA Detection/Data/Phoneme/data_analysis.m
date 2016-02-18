%Phoneme data analysis.
%The import_matlab.m script should be used to import the original text 
%phoneme.data file into phoneme.mat
load('phoneme.mat','phoneme','class');

%%
index = zeros(numel(class),1);
str = 'dcl'; %'aa';
for k = 1:numel(class)
    index(k) = strcmp(class{k},str);
end

subset = phoneme(index==1,:);
[n,p] = size(subset);
if strcmp(str,'dcl')
    dcl = subset;
    save('dcl.mat','dcl','n','p');
end
%%
%dcl: nxp, where n=757, p=256
%n instances of  phoneme "dcl" periodogram on 256 points.
type_arr = {'cov','cor'};
for i=1:length(type_arr)
    type = type_arr(i);
    type = type{:};
    switch type
        case 'cov'
            S = cov(subset);
        case 'cor'
            S = corr(subset);
    end
    lambda = eig(S);
    %% plot
    if strcmp(type,'cov')
        lambda = lambda/mean(lambda);
    end
    gamma = p/n;
    gamma_plus = (1+sqrt(gamma))^2;
    gamma_minus = (1-sqrt(gamma))^2;
    UB=10;
    num_large = sum(lambda> UB);
    %cov: num_large = 10, num_sig = 62
    %cov normalized: num_large = 2, num_sig = 9
    %corr: num_large = 2; num_sig = 6
    MP = @(x) sqrt(max(0,((1+sqrt(gamma))^2-x).*(x-(1-sqrt(gamma))^2)))./(2*pi*gamma.*x);
    grid = linspace(0,max(lambda(lambda<UB)),500);
    density = MP(grid);
    
    [heights,locations] = hist(lambda(lambda<UB),10*floor(sqrt(p)));
    width = locations(2)-locations(1);
    heights = heights / max(heights);
    
    %plot
    figure, hold on
    bar(locations,heights);
    plot(grid, density/max(density),'r','LineWidth',4)
    set(gca,'fontsize',20)
    xlabel('Eigenvalues');
    %h = legend( 'Empirical','MP','Location','best');
    %set(h,'FontSize',20);
    
    %%
    filename = sprintf( './Img/phoneme = %s type = %s n = %d p = %d.png',str,type, n,p);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
    %% How many significant?
    TW_UB = 4*n^(-1)* (sqrt(n-1)+sqrt(p))*(1/sqrt(n-1)+1/sqrt(p))^(1/3);
    num_sig = sum(lambda> (1+sqrt(gamma))^2+TW_UB);
end

%% Rescale sample cov
%%
type ='cov';
S = cov(subset);

lambda = eig(S);
sig_arr = linspace(0.5,5,6);
num_sig = zeros(length(sig_arr),1);
for i=1:length(sig_arr)
    sig = sig_arr(i);
    lambda = sig*lambda/mean(lambda);
    num_large = sum(lambda> UB);
    grid = linspace(0,max(lambda(lambda<UB)),500);
    density = MP(grid);
    
    [heights,locations] = hist(lambda(lambda<UB),10*floor(sqrt(p)));
    width = locations(2)-locations(1);
    heights = heights / max(heights);
    
    %plot
    figure, hold on
    bar(locations,heights);
    plot(grid, density/max(density),'r','LineWidth',4)
    set(gca,'fontsize',20)
    xlabel('Eigenvalues');
    
    %%
    filename = sprintf( './Img/Sigma scaling/phoneme = %s type = %s sig = %.2f n = %d p = %d.png',str,type,sig, n,p);
    saveas(gcf, filename,'png');
    fprintf(['Saved Results to ' filename '\n']);
    close(gcf)
    num_sig(i) = sum(lambda> (1+sqrt(gamma))^2+TW_UB);
    num_sig(i) = num_sig(i)+sum(lambda< (1-sqrt(gamma))^2-TW_UB);
end

%%
save('./Img/Sigma scaling/num_sig_dcl.mat','sig_arr','num_sig');

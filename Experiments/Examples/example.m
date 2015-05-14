%% examples for using the EigenEdge package, specifically, the Spectrode function
%% Quick example
t = [1; 5]; %location of population eigenvalues: H = 1/n\sum delta_{t_i}
gamma = 1/2; %aspect ratio gamma = p/n
[grid, density] =  spectrode(t, gamma); %compute limit spectrum
figure, plot(grid, density,'r','LineWidth',4) %plot
xlabel('Eigenvalue')
ylabel('Density');

%%
filename = sprintf( './quick_example.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
%% Long example:  mixture of uniform + arithmetic model

gap = 0.01;
t_min = 2;
t_max = 11;
K = 10;
[t,w]  = comb_model(K, t_min,t_max,gap);
w = 1/2*w;

r = [0.5 1.5];
w_int = 1/2;

epsi = 10^(-5);
gamma = 0.01;
[grid, density,~, ~, ~, K_hat, l_hat, u_hat] =  spectrode(t, gamma, w, r,w_int,epsi);

%empirical eigenvalues
rng(0);
n = 5*1e4;
p = floor(gamma*n);
pop_eigs = random_draw_disc_unif(t,w,r,w_int, p);
[X] = 1/sqrt(n)*randn(n,p)*diag(sqrt(pop_eigs));
D = svd(X,'econ').^2;
%
[heights,locations] = hist(D,10*floor(sqrt(p)));
width = locations(2)-locations(1);
heights = heights / max(heights);

%plot
figure, hold on
bar(locations,heights);
plot(grid, density/max(density),'r','LineWidth',4)
set(gca,'fontsize',14)
xlabel('Eigenvalues');
h = legend( 'Empirical Eigenvalues','Theoretical Prediction','Location','best');
set(h,'FontSize',14);
%%
filename = sprintf( './Illustration_mixture_3.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
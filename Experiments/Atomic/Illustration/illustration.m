%Test ODE method for solving MP equation

%% Mixture of three point masses
%theoretical lsd
t = (1:3)';
w = ones(length(t),1)/length(t);
gamma = 0.1;
[grid, density] =  compute_esd_ode(t,w,gamma);

%empirical eigenvalues
n = 1000;
p = floor(gamma*100);
testType = 'cluster';
[X] = generateCorrelatedData(n,p,testType);
D = svd(X,'econ').^2;
[heights,locations] = hist(D,3*floor(sqrt(p)));
width = locations(2)-locations(1);
heights = heights / max(heights);

%plot
figure
hold on
bar(locations,heights);
plot(grid, density/max(density),'r','LineWidth',4)
legend( 'Empirical Eigenvalues','Theoretical Prediction','Location','best');

%%
filename = sprintf( './Illustration_mixture.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
%%
figure;
w = ones(length(t),1)/length(t);
bar(t,w,0.5);
filename = sprintf( './Illustration_mixture_pop.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

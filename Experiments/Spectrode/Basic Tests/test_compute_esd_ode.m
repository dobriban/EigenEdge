%Test the ODE methods for solving MP equation
%this is a simple script with a few plots

%%
a = 0.5;
b = 1.5;
K = 6;
[t,w]  = geometric_model(a,b,K);
gamma = 2;
epsi  = 1e-6;
[grid,density,~,~,~,K_hat] =  compute_esd_ode(t,w,gamma,epsi,'brent');
%
figure, 
plot(grid, density/max(density),'r','LineWidth',4)
%% empirical eigenvalues
p = 100;
n = floor(p/gamma);
pop_eigs = random_discrete_draw(t,w,p);
[X] = 1/sqrt(n)*randn(n,p)*diag(sqrt(pop_eigs));
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

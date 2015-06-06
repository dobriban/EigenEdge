% Test fixed point iteration with complex starting point and real iterative
% update
addpath('C:\Code\EigenEdge\Code')
%%
gamma_plus = (1+sqrt(gamma))^2;
gamma_minus = (1-sqrt(gamma))^2;
MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;
%% Test correctness by computing theoretical and numerical solutions
gamma = 1/2;
t =1;
epsilon = 1e-4;
maxIter = 1e4;
grid = linspace(2,2,1)';
w = 1;
[density,m,v,numIter] = compute_esd_fp_test(t,w,gamma,epsilon,grid,maxIter);
theor_density = MP_density(grid);
err= 1/pi*imag(m)-theor_density;

%%
figure, hold on
ind = (grid>0.01)&(grid<1.1*(1+sqrt(gamma))^2*max(t));
plot(grid(ind),log10(abs(err(ind))),'linewidth',3,'color',rand(1,3))
xlabel('Eigenvalue')
ylabel('Log10 of error: theoretical vs numerical density');
figtitle('Iterative Method: Numerical error as a function of tolerance');

%% plot
%plotting the trajectory from compute_esd_fp_test
filename = sprintf( './real_iteration_trajectory.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
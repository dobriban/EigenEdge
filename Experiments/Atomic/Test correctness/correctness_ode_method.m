%Test ODE method for solving MP equation
%% Test correctness by computing theoretical and numerical solutions
% Null case where all eigenvalues are equal to 1
p = 200;
gamma = 1/2;
epsi_array = [4,6,8];
t = 1;
rng(0)
figure, hold on
a = {'--','-.','-','-.'};
for i=1:3;
    gamma_plus = (1+sqrt(gamma))^2;
    gamma_minus = (1-sqrt(gamma))^2;
    MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;
    epsilon = 10.^(-epsi_array(i));
    w = 1;
    [grid, density_ode,m, mass_at_0] = compute_esd_ode(t,w,gamma,epsilon);
    theor_density = MP_density(grid);
    err= density_ode-theor_density;
    
    %plot the error
    h = plot(grid,log10(abs(err)),'linewidth',3,'color',rand(1,3));
    set(h,'LineStyle',a{i});
end
set(gca,'fontsize',14)
h = legend('\epsilon=1e-4','\epsilon=1e-6','\epsilon=1e-8','location','Best');
xlabel('x_i')
ylabel('\Delta(x_i,\epsilon)');
%figtitle('Numerical error as a function of tolerance');
set(h,'FontSize',14);
%% plot
filename = sprintf( './Correctness_tests_numerical_error_null.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% Test correctness by computing theoretical and numerical solutions
% Mixture of two point masses
p = 200;
gamma = 1/10;
epsi_array = [4,6,8];
t = [1,8];
w = [1,1]/2;
tau = 8;
gamma_plus = (1+sqrt(gamma))^2;
gamma_minus = (1-sqrt(gamma))^2;
a = {'--','-.','-','-.'};
rng(0)
figure, hold on
for i=1:3;
    
    epsilon = 10.^(-epsi_array(i));
    [grid, density_ode, m, mass_at_0] = compute_esd_ode(t,w,gamma, epsilon);
    [m_theor] = SparseStieltjes2(grid,gamma,1/2,tau);
    theor_density = 1/pi*imag(m_theor);
    err= density_ode -theor_density;
    
    ind = (min(find(theor_density>0)):max(find((theor_density>0))));
    
    h = plot(grid(ind),log10(abs(err(ind))),'linewidth',3,'color',rand(1,3));
    set(h,'LineStyle',a{i});
end
set(gca,'fontsize',14)
h = legend('\epsilon=1e-4','\epsilon=1e-6','\epsilon=1e-8','location','Best');
xlabel('x_i')
ylabel('\Delta(x_i,\epsilon)');
%figtitle('Numerical error as a function of tolerance');
set(h,'FontSize',14);
%% plot
filename = sprintf( './Correctness_tests_numerical_error_two-point.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);
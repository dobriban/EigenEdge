%Test Iterative method for solving MP equation
%This code is very similar to the correctness tests for the ODE method
%% Test correctness by computing theoretical and numerical solutions
% Null case where all eigenvalues are equal to 1
p = 200;
gamma = 1/2;
epsi_array = linspace(1,2,3);
t =1;
rng(0)
figure, hold on
for i=1:length(epsi_array);
    gamma_plus = (1+sqrt(gamma))^2;
    gamma_minus = (1-sqrt(gamma))^2;
    MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;
    epsilon = 10.^(-epsi_array(i));
    
    maxIter = Inf;
    n = floor(p/gamma);
    grid = linspace(1e-3,10,1e3)';
    w = 1;
    tic
    [~,m,~] = mp_solve_iter(t,w,gamma,epsilon,grid,maxIter);
    t1 = toc;
    fprintf('Test %d/%d: For accuracy %2.2e, ITER method took %2.2f seconds\n',i,length(epsi_array), epsilon, t1);
   
    theor_density = MP_density(grid);
    err= 1/pi*imag(m)-theor_density;
    
    ind = (grid>0.01)&(grid<1.1*(1+sqrt(gamma))^2*max(t));
    
    %plot the two densities
    %     figure,
    %     plot(f(ind),1/pi*imag(m(ind))); hold on
    %     plot(f(ind),theor_density(ind),'rx');
    %
    %plot the error
    plot(grid(ind),log10(abs(err(ind))),'linewidth',3,'color',rand(1,3))
end
legend('epsi=1e-1','epsi=3e-2','epsi=1e-2','location','Best')
xlabel('Eigenvalue')
ylabel('Log10 of error: theoretical vs numerical density');
figtitle('Iterative Method: Numerical error as a function of tolerance');

%% plot
filename = sprintf( './Correctness_iter_tests_numerical_error_null.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);


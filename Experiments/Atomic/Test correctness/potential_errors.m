%an example where the ODE solver computes apparently
%erroneous solutions
%Reason: the accuracy parameter epsilon is not small enough 
t = [1;1;2;2;2;3;3;3;3;3];
p = 10;
n =200;
epsilon = 1e-3;
[f, m, mass_at_0] = mp_solve_ode(t,gamma,epsilon);
[~,m_iter] = mp_solve_iter(t,n,epsilon,f);

ind = (f>0.01)&(f<1.1*(1+sqrt(gamma))^2*max(t)); 
[heights,locations] = hist(t,3*floor(sqrt(p)));
width = locations(2)-locations(1);
heights = heights / max(heights);

figure
hold on
bar(locations,heights);
plot(f(ind), imag(m(ind))/max(imag(m(ind))),'r','LineWidth',2)
plot(f(ind), imag(m_iter(ind))/max(imag(m_iter(ind))),'g','LineWidth',2)
legend( 'Population Eigenvalues','ODE','Iterative', 'Location','best');

%% plot
filename = sprintf( './Potential_error_three-point-mixture.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% Try to solve the problem by setting smaller epsilon
epsilon = 1e-5;
[f, m] = mp_solve_ode(t,gamma,epsilon);
ind = (f>0.01)&(f<1.1*(1+sqrt(gamma))^2*max(t)); 

figure
hold on
bar(locations,heights);
plot(f(ind), imag(m(ind))/max(imag(m(ind))),'r','LineWidth',2)
plot(f(ind), imag(m_iter(ind))/max(imag(m_iter(ind))),'g','LineWidth',2)
legend( 'Population Eigenvalues','ODE','Iterative', 'Location','best');

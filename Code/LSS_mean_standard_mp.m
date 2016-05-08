function int = LSS_mean_standard_mp(g, gamma,epsilon)
%compute the mean of a LSS g in the standard MP model
%this is a contour integral 

%Used in the Examples of the Spectrode paper. 

if ~exist('epsilon','var')
    epsilon = 1e-10;
end

%% find starting point of ODE
a = 1.1*(1+sqrt(gamma))^2;

silv_fn = @(v) -1./v + gamma./(1+v)-a;

grid_neg = linspace(-1+1e2*sqrt(epsilon),-1e2*sqrt(epsilon),sqrt(floor(1/epsilon)));
z = silv_fn( grid_neg);
c_gamma = grid_neg(z==min(z));
c_gamma = c_gamma(1);

ode_start_point = fzero(silv_fn,[c_gamma,-epsilon]);

%plot(grid_neg,z)
%should be 0
%abs(silv_fn(ode_start_point))
%% solve ODE
b = a/2;
F = @(h) 1./(1./h.^2-gamma./(1+h).^2);
c = @(t) a/2 + b * exp(2*pi*1i*t);
dc =  @(t) 2*b*pi*1i* exp(2*pi*1i*t);
h_diff = @(t,h) F(h).*dc(t);
grid = linspace(0,1,sqrt(floor(1/epsilon)))';
ep = epsilon;
options = odeset('RelTol', ep, 'AbsTol',ep);

[~,h] = ode45(h_diff,grid,ode_start_point,options);

% figure, hold on
% plot(grid,real(h))
% plot(grid,imag(h))

%%
G = @(z,v) g(z).* gamma.* v.^3./(1+v).^3./[1 - gamma.* v.^2./(1+v).^2].^2;

c_arr=c(grid);
obj_fun = G(c_arr,h).*dc(grid);
int  = -1/(2*pi*1i)* trapz(grid,obj_fun);
int = real(int);




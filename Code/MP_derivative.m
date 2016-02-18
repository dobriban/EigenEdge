function [grid, density,ind_outside,ind_in,K_hat, l_hat, u_hat,a,b,upper_pt,v_grid] = MP_derivative(t,w_null,spikes,w_alt,gamma,epsi)
%compute the density of the weak derivative of the MP transform
%dF(H,G_0)
%where
%H = delta_1
%G_0 = G - delta_1
%with
%H = \sum_i  w_null(i) * delta(t(i))
%G = \sum_i  w_alt(i) * delta(spikes(i))

%Inputs
%t - null eigenvalues
%w_alt - null eigenvalue weights
%spikes - loc of spikes
%w_alt - spike weights
%gamma - aspect ratio
%grid - grid where to compute LSS
%
%Outputs
%grid - grid where LSS is computed
%density - density of the weak derivative
%ind_outside,ind_in - indices of the grid outside and inside the bulk
% K_hat, l_hat, u_hat - SPECTRODE outputs (support of ESD)
% a,b - intervals where v is increasing
% upper_pt - upper phase transition threshold

%epsi = 5*1e-4;
%the grid generated here has too large gaps at the edges
%the grid outside the bulk is super dense
%the grid inside the bulk is much coarser, so that at the edge we have: 
%dgrid=[1e-05; 5e-06; 0; 0.04; 0.04]
%the problem in  compute_esd_ode is that the 
%length of the support interval is too large compared to the number of grid
%points used within it
%in particular, I always use 
%M = floor(sqrt(1/epsilon))+3;
%grid points
%after the change:
%[1e-05; 5e-06; 0; 0.003; 0.003]
[~, ~, ~, ~, ~, K_hat, l_hat, u_hat,~,~,grid,v_grid,a,b] = compute_esd_ode(t, w_null, gamma,epsi);  
dgrid = grid(2:length(grid))-grid(1:length(grid)-1);
dgrid = [dgrid; dgrid(length(dgrid))];
ind_nonsing = (dgrid>0);
grid = grid(ind_nonsing);
v_grid = v_grid(ind_nonsing);
%get rid of large values in the grid
good_ind = (grid<1.1*(1+sqrt(gamma))^2*max(max(spikes),max(t)));
grid = grid(good_ind);
v_grid = v_grid(good_ind);

v_prime = @(v) (1/v^2 - gamma* sum( w_null.*(t.^(-1) + v).^(-2)))^(-1);
v_prime_grid =  zeros(length(v_grid),1);
for i=1:length(v_grid)
    v_prime_grid(i) = v_prime(v_grid(i));
end

numer = zeros(length(v_grid),1); 
for i=1:length(v_grid)
    numer(i) = sum(spikes.*w_alt./(1+spikes*v_grid(i)))- sum(t.*w_null./(1+t*v_grid(i)));
end
st = -v_prime_grid.*gamma.*numer;
density = 1/pi*imag(st);


%maybe I'm too conservative here?
ind_outside = ones(length(grid),1);
ind_outside_plus = ones(length(grid),1);
%epsi = 10*epsi;
L = 1;
%L = 0;
for i=1:K_hat
    ind_outside((grid>=l_hat(i))&(grid<=u_hat(i)))=0;
    %also set a small portion of the inside to 
    l = find(grid==l_hat(i));
    u = find(grid==u_hat(i));
    ind_outside_plus(min(l+L,u):max(u-L,l))=0;
end
ind_in = (ind_outside_plus==0);

%compute upper phase transition threshold
%in theory this equals x = -1/v(u)
%where u is upper edge of spectrum
u = max(u_hat);
v_u = interp1(grid, v_grid,u);
upper_pt = -1/v_u;


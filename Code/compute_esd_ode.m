function [grid, density, m, v, mass_at_0, K_hat, l_hat, u_hat, v_x, f_hat,grid_large,v_large,a,b] = ...
    compute_esd_ode(t, w, gamma,epsilon,edge_finding,M)
%Compute limit spectrum of covariance matrices with ODE method
%uses the inverse of the Stieltjes transform obtained by Silverstein&Choi
%(95) to determine the edges of the spectrum
%
%Then uses MP differential equation in an ODE solver
%
%Written by Edgar Dobriban
% Inputs
% t - population eigenvalues >=0
% gamma - aspect ratio p/n
% w - mixture weights >=0; default: w = uniform
% epsilon - square root of the imaginary part of the grid where MP is solved
% edge_finding - specify the method by which the edges are computed. 'grid'
% - default - computes the inverse ST on a grid; 'brent' - uses Brent's
% method to find increasing intervals of ST.
% M - number of grid points in each interval of the support

%Outputs
% grid - real grid where Stieltjes transform is evaluated, real vector of size Nx1
% density - approximation to the density on the grid
% m - numerical approximation to Stieltjes transform on real line,
%       complex vector of size Nx1
% v - numerical approximation to dual Stieltjes transform on real line,
%       complex vector of size Nx1
% mass_at_0 - discrete probability mass placed at 0 by the dual ST
% K_hat - numerical approximation to number of disjoint clusters of support
% l_hat - lower endpoints of support intervals; real vector of size K_hat;
% u_hat - upper endpoints of support intervals; real vector of size K_hat;
% v_x - grid points within support intervals where density is approximated;
%       K_hat x 1 cell array; cell i contains the i-th support interval
%        (these are a superset of the indices in grid; and they are
%        provided for convenience in some applications)
% f_hat - numerical approximation of density on grid v_x; same format as v_x
% grid_large,v_large - positive entries of the grid, with their v values
%   (these are supersets of grid and v; of the same length)
% a,b - intervals where v is increasing

if ~exist('epsilon','var')
    epsilon = 1e-4;
end
%set lower bound on epsilon
%the limitation comes from the ODE solver, which ends up working with 
%2x the precision of epsilon.
%and currently this is limited to 16 digits
epsilon = max(epsilon, 1e-8);

delta_v = sqrt(epsilon);

p = length(t);
if ~exist('w','var')
    w = 1/p*ones(p,1);
end

if ~exist('M','var')
    M = floor(sqrt(1/epsilon))+3;
end

if ~exist('edge_finding','var')
    edge_finding = 'grid';
    %edge_finding = 'brent';
end

options = optimset('TolX',10^(-20));
%compute mass of sample distribution at zero
p0 = sum(w.*(t==0));
mass_at_0 = max(0, 1-gamma^(-1)*(1-p0));

T = unique(t(t>0));
B = sort(-1./T);
num_unique = length(B);
z_l_endpoints = zeros(num_unique+2,1);
z_u_endpoints = zeros(num_unique+2,1);
z_double_prime_equiv = @(v) sum( w.*(t.*v).^3./(t.*v+1).^3)-1/gamma;
z_prime_equiv = @(v) sum( w.*(t.*v).^2./(t.*v+1).^2)-1/gamma;
z = @(v) - 1./v + gamma* sum(w.*t./(1 + t.*v));

%% outside the support
z_grid = zeros(num_unique+2,M); %grids in z-space
v_grid = zeros(num_unique+2,M); %grids in v-space
num_grid_points = zeros(num_unique+2,1); %number of points in each grid segment

%go interval by interval;
num_clus = 1;
z_l_endpoints(1) = 0;
if (gamma <1) %if need to look at lower half
    v_L = min(B) -1;
    %if epsilon is so large that the grid is empty, make it smaller
    while (v_L+epsilon)>=(min(B)-epsilon)
        epsilon = epsilon/2;
    end
    switch edge_finding
        case 'grid'
            v= linspace(v_L+epsilon, min(B)-epsilon, M);
            [grid,~,~,maxz,ind_max] = evaluate_inverse_ST(t,w,gamma, v);
            
            while (ind_max==find((grid<Inf)==1, 1 ))  %need to deal with Inf's
                v_L = 10*(v_L - min(B)) + min(B);
                v= linspace(v_L+epsilon, min(B)-epsilon, M);
                [grid,~,~,maxz,ind_max] = evaluate_inverse_ST(t,w,gamma, v);
            end
            z_u_endpoints(num_clus) = maxz; 
            num_grid_points(num_clus) = ind_max; %number of grid points in first interval
            z_grid(num_clus,1:ind_max) = grid(1:ind_max);
            v_grid(num_clus,1:ind_max) = v(1:ind_max);
        case 'brent'
            v_0 = [v_L+epsilon, min(B)-epsilon]; % initial interval
            while (z_prime_equiv(v_0(1))*z_prime_equiv(v_0(2))>0)
                v_L = 10*(v_L - min(B)) + min(B);
                v_0 = [v_L+epsilon, min(B)-epsilon]; % initial interval
            end
            v_x = fzero(z_prime_equiv,v_0,options); %v_L<v_x<v_U, and [v_L,v_x] is interesting region
            z_u_endpoints(num_clus) = z(v_x);
            v_local = linspace(v_L,v_x,M);
            [num_grid_points(num_clus), z_grid(num_clus,1:M), v_grid(num_clus,1:M)] = grid_output(v_local,z);
    end
else
    num_grid_points(num_clus) = 0;
end

%each interval between point masses -1/t
for i = 1:num_unique-1
    v_L = B(i); %lower endpoint -1/t
    v_U = B(i+1); %upper endpoint -1/t
    %if epsilon is so large that the grid is empty, make it smaller
    while (v_L+epsilon)>=(v_U-epsilon)
        epsilon = epsilon/2;
    end
    %Find the edges of the spectrum
    switch edge_finding
        case 'grid'
            v= linspace(v_L+epsilon, v_U-epsilon, M); %adaptive grid-forming
            [grid, ~, ~,~,~,decreasing,local_min,local_max, local_min_ind, local_max_ind] ...
                =  evaluate_inverse_ST(t,w,gamma, v);
            if (decreasing==0)
                num_clus = num_clus + 1;
                z_l_endpoints(num_clus) = local_min;
                z_u_endpoints(num_clus) = local_max;
                num_grid_points(num_clus) = local_max_ind - local_min_ind + 1; %num grid points where increasing
                c = num_grid_points(num_clus);
                z_grid(num_clus,1:c) = grid(local_min_ind:local_max_ind);
                v_grid(num_clus,1:c) = v(local_min_ind:local_max_ind);
            end
        case 'brent'
            v_0 = [v_L+delta_v, v_U-delta_v]; % initial interval
            while (z_double_prime_equiv(v_0(1))*z_double_prime_equiv(v_0(2))>0)
                delta_v = delta_v/2;
                v_0 = [v_L+delta_v, v_U-delta_v];
            end
            m0 = max(z_double_prime_equiv(v_0(1)), z_double_prime_equiv(v_0(2)));
            z_2_n = @(v) z_double_prime_equiv(v)/m0;
            v_x = fzero(z_2_n,v_0,options);
            if (z_prime_equiv(v_x)<0)
                num_clus = num_clus + 1;
                v_0 = [v_L+delta_v v_x]; %to the left of the point where z''=0
                while (z_prime_equiv(v_0(1))*z_prime_equiv(v_0(2))>0)
                    delta_v = delta_v/2;
                    v_0 =  [v_L+delta_v v_x];
                end
                left_zero = fzero(z_prime_equiv,v_0,options);
                z_l_endpoints(num_clus) = z(left_zero);
                v_0 = [v_x v_U-delta_v]; %to the right of the point where z''=0
                while (z_prime_equiv(v_0(1))*z_prime_equiv(v_0(2))>0)
                    delta_v = delta_v/2;
                    v_0 = [v_x v_U-delta_v];
                end
                right_zero = fzero(z_prime_equiv,v_0,options);
                z_u_endpoints(num_clus) = z(right_zero);
                v_local = linspace(left_zero,right_zero,M);  %possible bug here
                %v_local = (z_l_endpoints(num_clus):v_delta:z_u_endpoints(num_clus))';
                [num_grid_points(num_clus), z_grid(num_clus,1:M), v_grid(num_clus,1:M)] = grid_output(v_local,z);
            end
    end
end

%last interval between -1/t & 0
v_L = B(num_unique);
v_U = 0;
%if epsilon is so large that the grid is empty, make it smaller
while (v_L+epsilon)>=(v_U-epsilon)
    epsilon = epsilon/2;
end
switch edge_finding
    case 'grid'
        v= linspace(v_L+epsilon, v_U-epsilon, M); %adaptive grid-forming
        [grid, minf,ind_min] = evaluate_inverse_ST(t,w,gamma, v);
        num_clus = num_clus + 1;
        z_l_endpoints(num_clus)  = minf;
        z_u_endpoints(num_clus)  = Inf;
        c = M-ind_min+1;
        num_grid_points(num_clus) = c; %num grid points where increasing
        z_grid(num_clus,1:c) = grid(ind_min:M);
        v_grid(num_clus,1:c) = v(ind_min:M);
        
    case 'brent'
        v_0 = [v_L+epsilon v_U-epsilon]; % initial interval in v
        while (z_prime_equiv(v_0(1))*z_prime_equiv(v_0(2))>0)
            v_L = 10*(v_L - min(B)) + min(B);
            v_0 = [v_L+epsilon, min(B)-epsilon]; % initial interval
        end
        v_x = fzero(z_prime_equiv,v_0,options); %value of v where z'(v0)=0
        num_clus = num_clus + 1;
        z_l_endpoints(num_clus) = z(v_x); %z(v0)
        z_u_endpoints(num_clus) = Inf;
        v_local = linspace(v_x,v_U,M)';
        [num_grid_points(num_clus), z_grid(num_clus,1:M), v_grid(num_clus,1:M)] = grid_output(v_local,z);
end

%the positive line: between 0 and inf
num_clus = num_clus+1;
z_l_endpoints(num_clus) = -Inf;
if (gamma >1) %if upper half maximum may exceed 0
    v_U = 1;
    %if epsilon is so large that the grid is empty, make it smaller
    while (epsilon)>=(v_U-epsilon)
        epsilon = epsilon/2;
    end
    switch edge_finding
        case 'grid'
            v= linspace(epsilon, v_U-epsilon, M);
            [grid,~,~,maxz,ind_max] = evaluate_inverse_ST(t,w,gamma, v);
            while (ind_max==length(v))  %no need to deal with Inf's
                v_U = 10*v_U;
                v= linspace(epsilon, v_U-epsilon, M);
                [grid,~,~,maxz,ind_max] = evaluate_inverse_ST(t,w,gamma, v);
            end
            z_u_endpoints(num_clus) =  maxz;
            num_grid_points(num_clus) = ind_max; %number of grid points in last interval
            z_grid(num_clus,1:ind_max) = grid(1:ind_max);
            v_grid(num_clus,1:ind_max) = v(1:ind_max);
            
        case 'brent'
            v_0 = [epsilon, v_U-epsilon];
            while (z_prime_equiv(v_0(1))*z_prime_equiv(v_0(2))>0)
                v_U = 10*v_U;
                v_0 = [epsilon, v_U-epsilon];
            end
            v_x = fzero(z_prime_equiv,v_0,options);
            z_u_endpoints(num_clus) = z(v_x);
            v_local = linspace(epsilon,v_x,M)';
            [num_grid_points(num_clus), z_grid(num_clus,1:M), v_grid(num_clus,1:M)] = grid_output(v_local,z);
    end
else %if gamma<1, this interval gives me the whole region m>0;
    %know that the function is strictly increasing in this case
    v_U = M;
    %if epsilon is so large that the grid is empty, make it smaller
    while (epsilon)>=(v_U-epsilon)
        epsilon = epsilon/2;
    end
    v= linspace(epsilon, v_U-epsilon, M);
    [grid,~,~,~,~] = evaluate_inverse_ST(t,w,gamma, v);
    z_u_endpoints(num_clus) = 0;
    num_grid_points(num_clus) =M; %number of grid points in last interval
    z_grid(num_clus,1:M) = grid(1:M);
    v_grid(num_clus,1:M) = v(1:M);
end

%finish up by retaining only the nonzero clusters
z_l_endpoints = z_l_endpoints(1:num_clus);
z_u_endpoints = z_u_endpoints(1:num_clus);
num_grid_points =  num_grid_points(1:num_clus);
z_grid = z_grid(1:num_clus,:);
v_grid = v_grid(1:num_clus,:);

%% inside the support
%Now stitch together all intervals where Stieltjes transform is real, and
%fill in the intervals where the Stieltjes transform is complex-valued
% First Stitch together the real ST
% pre-allocate grids
% imag points + real points

%changed: do not preallocate grid---because we are forming it adaptively
%grid_length = (num_clus-2)*M + sum(num_grid_points);
%v = zeros(grid_length,1);
%grid = zeros(grid_length,1);
K_hat = num_clus-2;
l_hat =zeros(K_hat,1);
u_hat =zeros(K_hat,1);
%v_x = zeros(M,K_hat);
%f_hat = zeros(M,K_hat);
v_x = cell(K_hat,1);
f_hat =cell(K_hat,1);
a = zeros(K_hat+1,1);
b = zeros(K_hat+1,1);

grid(1:num_grid_points(num_clus)) = z_grid(num_clus, 1:num_grid_points(num_clus));
v(1:num_grid_points(num_clus)) = v_grid(num_clus, 1:num_grid_points(num_clus));
%store the increasing intervals
ind_ep =1;
a(ind_ep) = v_grid(num_clus, 1);
%b(ind_ep) = v_grid(num_clus, num_grid_points(num_clus));
b(ind_ep) = Inf;
ind = num_grid_points(num_clus); %current index

%if gamma>1 need to set the u-endpoint to higher than 0
if (gamma>1)
    z_u_endpoints(1) = max(grid(1:num_grid_points(num_clus)));
end

%if gamma <1 need to add in the first interval
if (gamma < 1)
    i = 1;
    grid(ind + 1:ind+num_grid_points(i)) = z_grid(i, 1:num_grid_points(i));
    v(ind + 1:ind+num_grid_points(i)) = v_grid(i, 1:num_grid_points(i));
    ind_ep  = ind_ep+1;
    a(ind_ep) = -Inf;
    b(ind_ep) = v_grid(i, num_grid_points(i));
    ind = ind + num_grid_points(i);
end

%store the increasing intervals
for i=2:num_clus-1
    ind_ep =ind_ep+1;
    a(ind_ep) = v_grid(i, 1);
    b(ind_ep) = v_grid(i, num_grid_points(i));
end

%functions for ODE solving
MP_diff = @(time,v) (1/v^2 - gamma* sum( w.*(t.^(-1) + v).^(-2)))^(-1);
ep = epsilon;
options = odeset('RelTol', ep, 'AbsTol',ep);

%from now on, can handle all points in a uniform way
%have num_clus-2 pairs of intervals left: look at the support first, then
%outside
for i=2:num_clus-1
    %within the support
    %define the parameters of an interval in the support
    endpoint_1 = z_u_endpoints(i-1); %lower endpoint in f-space
    endpoint_2 = z_l_endpoints(i); %upper endpoint in f-space
    l_hat(i-1) = endpoint_1;
    u_hat(i-1) = endpoint_2;
    
    %the grid within the i-th support interval
    %adaptively set the grid length to have an accuracy approximately
    %sqrt(ep) within the support interval
    M_curr = floor((endpoint_2-endpoint_1)/sqrt(epsilon))+1;
    grid_current = linspace(endpoint_1,endpoint_2,M_curr)';
    v_x(i-1) = {grid_current};
    
    %find starting point for ODE using iterative method:
    [~,~,v_start] =  compute_esd_fp(t,w,gamma,epsilon, grid_current(1));
    
    %find dual Stieltjes transform using ode
    [~,v0] = ode45(MP_diff,grid_current,v_start,options);
    
    %set the output
    grid(ind+1:ind + M_curr) = grid_current;
    v(ind+1:ind + M_curr) = v0;
    f_hat(i-1) = {1/pi*imag( 1/gamma*v0- (1-1/gamma)./grid_current)};
    ind = ind + M_curr;
    
    %add in the corresponding interval outside the suppost
    grid(ind + 1:ind+num_grid_points(i)) = z_grid(i, 1:num_grid_points(i));
    v(ind + 1:ind+num_grid_points(i)) = v_grid(i, 1:num_grid_points(i));
    ind = ind + num_grid_points(i);
end

%retain the relevant part of the grid
good_ind = (grid>0)&(grid<(1+sqrt(gamma))^2*max(t)*1.1);
grid_large = grid(grid>0);
v_large = v(grid>0);
grid = grid(good_ind);
v = v(good_ind);
v=transpose(v);
v_large=transpose(v_large);

%remove negative jumps in the grid
ind_neg = ((grid(2:length(grid))-grid(1:length(grid)-1))<=0);
grid = grid(~ind_neg);
v = v(~ind_neg);

%compute the usual ST from the companion ST
m = 1/gamma*v- (1-1/gamma)./grid;
density = abs(1/pi*imag(m));


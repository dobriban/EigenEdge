function [grid, density, m, v, mass_at_0, K_hat, l_hat, u_hat, x, f_hat] = ...
    compute_esd_ode_non_atomic(t, w, gamma, r,w_int,epsilon,M)
%Compute limit spectrum of covariance matrices with ODE method
%uses the inverse of the Stieltjes transform obtained by Silverstein&Choi
%(95) to determine the edges of the spectrum
%Then uses MP differential equation in an ODE solver

% The input spectrum is a mixture
% H = sum w(i)*delta_t(i) + sum w_int(i)*unif[r(i,1), r(i,2)]
% where
% delta_a  - point mass at a
% unif(a,b) - uniform distribution on [a,b]

%Written by Edgar Dobriban
% Inputs
% t - population eigenvalues >0; real vector of size J,
% w - mixture weights >0; real vector of size J, default: w = uniform
% gamma - aspect ratio p/n, real scalar
% epsilon - (optional) accuracy parameter
% r - (optional) intervals where the continuous part of the spectrum is supported,
%       real matrix of size Jx2, i-th row encodes the endpoints [a,b] of
%       the i-th continuous component of the spectrum
% w_int - (optional) weights of the continuous part of the spectrum
% M - (optional) number of grid points in each interval of the support

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
% x - grid points within support intervals where density is approximated;
%       real matrix of size M x K_hat; column i contains the i-th support
%       interval
% f_hat - numerical approximation of density on grid x; same format as x

if ~exist('epsilon','var')
    epsilon = 1e-4;
end

p = length(t);
if ~exist('w','var')
    w = 1/p*ones(p,1);
end

if ~exist('M','var')
    M = floor(sqrt(1/epsilon))+3;
end

%is there a continuous component to the spectrum?
if ~exist('r','var')
    r = []; w_int = [];
end


%compute mass of sample distribution at zero
mass_at_0 = max(0, 1-gamma^(-1));

%find unique elements in the point masses and endpoints
T = unique([t; r(:,1); r(:,2)]);
B = sort(-1./T);
num_unique = length(B);

%exclude intervals covered by the continuous part of the support
excluded = zeros(num_unique,1);
for i=1:size(r,1)
    a = r(i,1); b = r(i,2);
    if (a>=b)
        fprintf('Continuous support intervals must have a<b\n');
    end
    ind_min=find(B>=-1/a,1,'first');
    ind_max=find(B<=-1/b,1,'last');
    excluded(ind_min:ind_max-1)=1;
end


l_endpoints = zeros(num_unique+2,1);
u_endpoints = zeros(num_unique+2,1);

%% outside the support
f_grid = zeros(num_unique+2,M); %grids in f-space
m_grid = zeros(num_unique+2,M); %grids in m-space
num_grid_points = zeros(num_unique+2,1); %number of points in each grid segment

%go interval by interval;
num_clus = 1;
l_endpoints(1) = 0;
if (gamma <1) %if need to look at lower half
    L = min(B) -1;
    v= linspace(L+epsilon, min(B)-epsilon, M);
    [f,~,~,maxf,ind_max] = evaluate_inverse_ST_non_atomic(t,w,r,w_int, gamma, v);
    
    while (ind_max==find((f<Inf)==1, 1 ))  %need to deal with Inf's
        L = 10*(L - min(B)) + min(B);
        v= linspace(L+epsilon, min(B)-epsilon, M);
        [f,~,~,maxf,ind_max] = evaluate_inverse_ST_non_atomic(t,w,r,w_int, gamma, v);
    end
    %plot(m,f) % should see a function with a unique maximum <0
    u_endpoints(num_clus) = maxf;
    num_grid_points(num_clus) = ind_max; %number of grid points in first interval
    f_grid(num_clus,1:ind_max) = f(1:ind_max);
    if ~issorted(f_grid(num_clus,1:ind_max))
        fprintf('Unsorted f grid in lower half.\n');
    end
    m_grid(num_clus,1:ind_max) = v(1:ind_max);
else
    %number of grid points in first interval = 0; worry not about setting
    %f_grid, m_grid
    num_grid_points(num_clus) = 0;
    %but what is u_endpoints(1)?
    %u_endpoints(1) = ? should be 0 i guess
end

%each interval between point masses -1/t
for i = 1:num_unique-1
    %only check it if it's not excluded
    if (excluded(i)==0)
        L = B(i);
        U = B(i+1);
        delta = (U-L)/(M+2);
        v= linspace(L+delta, U-delta, M); %adaptive grid-forming
        [f, ~, ~,~,~,decreasing,local_min,local_max, local_min_ind, local_max_ind] ...
            =  evaluate_inverse_ST_non_atomic(t,w,r,w_int, gamma, v);
        %plot(m,f)
        if (decreasing==0)
            num_clus = num_clus + 1;
            l_endpoints(num_clus) = local_min;
            u_endpoints(num_clus) = local_max;
            num_grid_points(num_clus) = local_max_ind - local_min_ind + 1; %num grid points where increasing
            c = num_grid_points(num_clus);
            f_grid(num_clus,1:c) = f(local_min_ind:local_max_ind);
            if ~issorted(f_grid(num_clus,1:c))
                fprintf('Unsorted f grid between clusters.\n');
            end
            m_grid(num_clus,1:c) = v(local_min_ind:local_max_ind);
        end
    end
end

%last interval between -1/t & 0
L = B(num_unique);
U = 0;
delta = (U-L)/(M+2);
v= linspace(L+delta, U-delta, M); %adaptive grid-forming
[f, minf,ind_min] = evaluate_inverse_ST_non_atomic(t,w,r,w_int, gamma, v);
num_clus = num_clus + 1;
l_endpoints(num_clus)  = minf;
u_endpoints(num_clus)  = Inf;
c = M-ind_min+1;
c = c(1);
num_grid_points(num_clus) = c; %num grid points where increasing
f_grid(num_clus,1:c) = f(ind_min:M);
if ~issorted(f_grid(num_clus,1:c))
    fprintf('Unsorted f grid in interval between -1/t & 0.\n');
end
m_grid(num_clus,1:c) = v(ind_min:M);
%plot(m,f)


%the positive line: between 0 and inf
num_clus = num_clus+1;
l_endpoints(num_clus) = -Inf;
if (gamma >1) %if upper half maximum may exceed 0
    U = 1;
    v= linspace(epsilon, U-epsilon, M);
    [f,~,~,maxf,ind_max] = evaluate_inverse_ST_non_atomic(t,w,r,w_int, gamma, v);
    
    while (ind_max==length(v))  %no need to deal with Inf's
        U = 10*U;
        v= linspace(epsilon, U-epsilon, M);
        [f,~,~,maxf,ind_max] = evaluate_inverse_ST_non_atomic(t,w,r,w_int, gamma, v);
    end
    %plot(m,f) % should see a function with a unique maximum > 0
    u_endpoints(num_clus) = maxf;
    num_grid_points(num_clus) = ind_max; %number of grid points in last interval
    f_grid(num_clus,1:ind_max) = f(1:ind_max);
    if ~issorted(f_grid(num_clus,1:ind_max))
        fprintf('Unsorted f grid in interval 0 & Inf, gamma>1.\n');
    end
    m_grid(num_clus,1:ind_max) = v(1:ind_max);
else %if gamma<1, this interval gives me the whole region m>0;
    %know that the function is strictly increasing in this case
    U = M;
    v= linspace(epsilon, U-epsilon, M);
    [f,~,~,~,~] = evaluate_inverse_ST_non_atomic(t,w,r,w_int, gamma, v);
    u_endpoints(num_clus) = 0;
    num_grid_points(num_clus) =M; %number of grid points in last interval
    f_grid(num_clus,1:M) = f(1:M);
    if ~issorted(f_grid(num_clus,1:M))
        fprintf('Unsorted f grid in interval 0 & Inf, gamma <1.\n');
    end
    m_grid(num_clus,1:M) = v(1:M);
    %plot(m,f) % should see a function with strictly increasing
end

%finish up by retaining only the nonzero clusters
l_endpoints = l_endpoints(1:num_clus);
u_endpoints = u_endpoints(1:num_clus);
num_grid_points =  num_grid_points(1:num_clus);
f_grid = f_grid(1:num_clus,:);
m_grid = m_grid(1:num_clus,:);

%% inside the support
%Now stitch together all intervals where Stieltjes transform is real, and
%fill in the intervals where the Stieltjes transform is complex-valued
% First Stitch together the real ST
% pre-allocate grids
% imag points + real points
grid_length = (num_clus-2)*M + sum(num_grid_points);
v = zeros(grid_length,1);
grid = zeros(grid_length,1);
K_hat = num_clus-2;
l_hat =zeros(K_hat,1);
u_hat =zeros(K_hat,1);
x = zeros(M,K_hat);
f_hat = zeros(M,K_hat);

grid(1:num_grid_points(num_clus)) = f_grid(num_clus, 1:num_grid_points(num_clus));
v(1:num_grid_points(num_clus)) = m_grid(num_clus, 1:num_grid_points(num_clus));
ind = num_grid_points(num_clus); %current index

%if gamma>1 need to set the u-endpoint to higher than 0
if (gamma>1)
    u_endpoints(1) = max(f(1:num_grid_points(num_clus)));
end

%if gamma <1 need to add in the first interval
if (gamma < 1)
    i = 1;
    grid(ind + 1:ind+num_grid_points(i)) = f_grid(i, 1:num_grid_points(i));
    v(ind + 1:ind+num_grid_points(i)) = m_grid(i, 1:num_grid_points(i));
    ind = ind + num_grid_points(i);
end

%functions for ODE solving
a = []; b = [];
if ~isempty(r)
    b = r(:,2);
    a = r(:,1);
end
conts_fun_diff = @(v) 1/v^2 - 1./(b-a).*(2*log((b.*v+1)./(a.*v+1)) + 1./(b.*v+1) - 1./(a.*v+1))/v^3;

MP_diff = @(time,v) 1/(1/v^2 - gamma*( sum( w.*(t.^(-1) + v).^(-2)) + sum(w_int.*conts_fun_diff(v)) ));
ep = epsilon;
options = odeset('RelTol', ep, 'AbsTol',ep);

%from now on, can handle all points in a uniform way
%have num_clus-2 pairs of intervals left: look at the support first, then
%outside
for i=2:num_clus-1
    %within the support
    %define the parameters of an interval in the support
    endpoint_1 = u_endpoints(i-1); %lower endpoint in f-space
    endpoint_2 = l_endpoints(i); %upper endpoint in f-space
    l_hat(i-1) = endpoint_1;
    u_hat(i-1) = endpoint_2;
    
    %the grid within the i-th support interval
    grid_current = linspace(endpoint_1,endpoint_2,M)';
    x(:,i-1) = grid_current;
    
    %find starting point for ODE using iterative method:
    [~,~,v_start] =  compute_esd_fp_non_atomic(t,w,r,w_int,gamma,epsilon, grid_current(1));
    
    %find dual Stieltjes transform using ode
    [~,v0] = ode45(MP_diff,grid_current,v_start,options);
    
    %set the output
    grid(ind+1:ind + M) = grid_current;
    v(ind+1:ind + M) = v0;
    f_hat(:,i-1) = 1/pi*imag( 1/gamma*v0- (1-1/gamma)./grid_current);
    ind = ind + M;
    
    %add in the corresponding interval outside the suppost
    grid(ind + 1:ind+num_grid_points(i)) = f_grid(i, 1:num_grid_points(i));
    v(ind + 1:ind+num_grid_points(i)) = m_grid(i, 1:num_grid_points(i));
    ind = ind + num_grid_points(i);
end

%retain the relevant part of the grid
maxi = max(max(t),max(r(:,2)));
good_ind = (grid>0)&(grid<(1+sqrt(gamma))^2*maxi*1.1);
grid = grid(good_ind);
v = v(good_ind);

%compute the usual ST from the dual ST
m = 1/gamma*v- (1-1/gamma)./grid;
density = abs(1/pi*imag(m));


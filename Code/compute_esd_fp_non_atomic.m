function [density,m,v,numIter,lastStepSize,v_d] = ...
    compute_esd_fp_non_atomic(t,w,r,w_int,gamma,epsilon, grid,multiplier_num_iter, maxIter,starting_point)
%Compute limit spectrum of covariance matrices with Fixed Point  method

%Inputs
%t - population eigenvalues, non-negative vector of length p
% w - mixture weights >=0; default: w = uniform
% r - intervals where the continuous part of the spectrum is supported,
%       real matrix of size Jx2, i-th row encodes the endpoints [a,b] of
%       the i-th continuous component of the spectrum
% w_int - weights of the continuous part of the spectrum
%gamma - aspect ratio p/n
%epsilon - size of imaginary part of grid
%grid - real grid where MP transform should be computed

%outputs
% grid - real grid where Stieltjes transform is evaluated, real vector of size Nx1
% density - approximation to the density on the grid
% m - numerical approximation to Stieltjes transform on real line, complex vector of size Nx1
% v - numerical approximation to dual Stieltjes transform on real line, complex vector of size Nx1
%v_d - estimate of the derivatves of the Gram Stieltjes transform of f on the grid
%numIter - number of iterations taken by the algorithm for each element of
%   the grid
%stepSize - last stepsize taken by the algorithm for each element of
%   the grid

%note that the tolerance is set to delta = 1e-8 in 'MP_transform'
%I know that these parameters guarantee that the result obeys
% |hat(v) - v| = O( epsilon^{1/2} + delta)
%thus epsilon is "more valuable"
%in the range of parameters of epsi ~ 1e-3; epsi is clearly the bottleneck

if ~exist('starting_point','var')
    starting_point = 'warm_start';
end
if ~exist('multiplier_num_iter','var')
    multiplier_num_iter = 1e2;
end
if ~exist('maxIter','var')
    maxIter = multiplier_num_iter/epsilon;
end
%set the convergence tolerance
if ~exist('tol','var')
    tol = epsilon;
end


grid_imag = grid + 1i*epsilon^2;
L = length(grid);
v = zeros(L,1);
v_d = zeros(L,6);
numIter = zeros(L,1);
lastStepSize = zeros(L,1);

a = []; b = [];
if ~isempty(r)
    b = r(:,2);
    a = r(:,1);
end

conts_fun = @(v) 1/v - 1./(b-a).*log((b.*v+1)./(a.*v+1))/v^2;
for i=1:L
    %define function
    z = grid_imag(i);
    fun = @(v) - 1/(z-gamma*(sum( w.*t./(1+t*v))+sum(w_int.*conts_fun(v))) );
    
    %starting point
    if (i==1)
        v1 = -1/grid_imag(i);
    else
        switch starting_point
            case 'default'
                v1 = -1/grid_imag(i);
            case 'warm_start'
                v1 = v(i-1);
        end
    end
    v1 = [v1 fun(v1(1))];
    j=2;
    
    while (abs(v1(j)-v1(j-1))>tol)&&(j<maxIter)
        v1 = [v1 fun(v1(j))];
        j = j+1;
    end
    v(i) = v1(j);
    numIter(i) = j;
    lastStepSize(i) = abs(v1(j)-v1(j-1));
end
m = 1/gamma*v- (1-1/gamma)./grid_imag;
density= 1/pi*imag(m);


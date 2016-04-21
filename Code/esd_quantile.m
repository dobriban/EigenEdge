function y = esd_quantile(t, w, gamma, q, epsilon)
%Compute a quantile of the limit spectrum of covariance matrices

% Written by Edgar Dobriban
% Inputs
% t - population eigenvalues >0
% w - mixture weights >0; default: w = uniform
% gamma - aspect ratio p/n
% q - quantile in (0,1) to compute
% epsilon - (optional) accuracy parameter

%Outputs
% y - approximate quantile of ESD
if ~exist('epsilon','var')
    epsilon = 1e-4;
end

[grid,density, ~, ~, mass_at_0] = compute_esd_ode(t, w, gamma,epsilon);

if length(q)==1
    y = esd_quantile_grid(grid,density,mass_at_0,q);
else
    y = esd_all_quantiles_grid(grid,density,mass_at_0,q);
end
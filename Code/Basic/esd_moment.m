function y = esd_moment(t, w, gamma, fun, epsilon)
%Compute a moment of the limit spectrum of covariance matrices

% Written by Edgar Dobriban
% Inputs
% t - population eigenvalues >0
% w - mixture weights >0; default: w = uniform
% gamma - aspect ratio p/n
% fun - function handle for moment to compute e.g. f = @(x) x; for mean
% epsilon - (optional) accuracy parameter

%Outputs
% y - approximate moment of ESD
if ~exist('epsilon','var')
    epsilon = 1e-4;
end

[grid,density, ~, ~, mass_at_0] = compute_esd_ode(t, w, gamma,epsilon);

y = esd_moment_grid(grid,density,mass_at_0,fun);

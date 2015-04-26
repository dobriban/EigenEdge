function y = esd_mode(t, w, gamma, epsilon)
%Compute the mode (or modes) of the limit spectrum of covariance matrices

% Written by Edgar Dobriban
% Inputs
% t - population eigenvalues >0
% w - mixture weights >0; default: w = uniform
% gamma - aspect ratio p/n
% epsilon - (optional) accuracy parameter

%Outputs
% y - approximate modes of ESD, a vector of real values
if ~exist('epsilon','var')
    epsilon = 1e-4;
end

[grid,density] = compute_esd_ode(t, w, gamma,epsilon);

y = esd_mode_grid(grid,density);
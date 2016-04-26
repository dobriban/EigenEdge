function y = esd_moment(t, w, gamma, fun, epsilon,companion)
%Compute a moment of the limit spectrum of covariance matrices

% Written by Edgar Dobriban
% Inputs
% t - population eigenvalues >0
% w - mixture weights >0; default: w = uniform
% gamma - aspect ratio p/n
% fun - function handle for moment to compute e.g. f = @(x) x; for mean
% epsilon - (optional) accuracy parameter
% companion - integrate with repect to companion (1) or usual (0) MP law
% (default=0 , usual MP)

%Outputs
% y - approximate moment of ESD
if ~exist('epsilon','var')
    epsilon = 1e-4;
end
if ~exist('companion','var')
    companion = 0;
end

[grid,density, ~, v, mass_at_0] = compute_esd_ode(t, w, gamma,epsilon);

if companion==0
    y = esd_moment_grid(grid,density,mass_at_0,fun);
else
    %rescale everything to companion
    density_companion = 1/pi*imag(v);
    if gamma>1
        mass_at_0_companion = 0;
    else
        %mass_at_0 = max(0, 1-gamma^(-1)*(1-p0));
        p0 = sum(w.*(t==0));
        mass_at_0_companion = max(0, 1-gamma*(1-p0));
    end 
    y = esd_moment_grid(grid,density_companion,mass_at_0_companion,fun);
end


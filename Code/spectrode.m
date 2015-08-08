function  [grid, density,m, v, mass_at_0, K_hat, l_hat, u_hat, x, f_hat]=...
    spectrode(t,gamma,varargin)
% Main interface to the Spectrode function for computing limit spectra of covariance matrices

% Example Usage:
% [grid, density] = Spectrode(t,gamma)
% where:
% t - a vector population eigenvalues >=0
% gamma - aspect ratio p/n

% For detailed usage, see attached documentation.
% [grid, density] = Spectrodes(t,gamma, `parameterName', `parameterValue',...) specifies additional parameters and options.


% The input spectrum is a mixture
% H = sum w(i)*delta_t(i) + sum w_int(i)*unif[r(i,1), r(i,2)]
% where
% delta_a  - point mass at a
% unif(a,b) - uniform distribution on [a,b]

% Inputs
% t - population eigenvalues >0; real vector of size J,
% gamma - aspect ratio p/n, real scalar
% w - mixture weights >0; real vector of size J, default: w = uniform
% r -  intervals where the continuous part of the spectrum is supported,
%       real matrix of size Jx2, i-th row encodes the endpoints [a,b] of
%       the i-th continuous component of the spectrum
% w_int -  weights of the continuous part of the spectrum
% epsilon - accuracy parameter
% M - number of grid points in each interval of the support

%Outputs
% grid - real grid where Stieltjes transform is evaluated, real vector of size Nx1
% density - approximation to the density on the grid
% m - numerical approximation to Stieltjes transform on real line,
%       complex vector of size Nx1
% v - numerical approximation to dual Stieltjes transform on real line,
%       complex vector of size Nx1

% Outputs only for Spectrode:
% mass_at_0 - discrete probability mass placed at 0 by the ESD
% K_hat - numerical approximation to number of disjoint clusters of support
% l_hat - lower endpoints of support intervals; real vector of size K_hat;
% u_hat - upper endpoints of support intervals; real vector of size K_hat;
% x - grid points within support intervals where density is approximated;
%       real matrix of size M x K_hat; column i contains the i-th support
%       interval
% f_hat - numerical approximation of density on grid x; same format as x

%% Process inputs
pa = parseArgumentsSpectrode(varargin{:});

parse(pa,t,gamma,varargin{:});
q = pa.Results;

%% Set Defaults
p = length(q.t);
if isempty(q.w)
    if ~isempty(q.r)
          error('The weights w must be specified if the continuous component is nonzero\n','missing_weights');
    else
    q.w = 1/p*ones(p,1);
    end
end

p = length(q.r);
if isempty(q.w_int)
    q.w_int = (1-sum(q.w))/p*ones(p,1);
end
%set default return arguments
K_hat = NaN; l_hat = NaN; u_hat = NaN; x = NaN; f_hat = NaN;

%if q.verbose <=1
%    warning('off','all');
%end

%% Perform computation

switch q.alg
    case 'spectrode'
        if isempty(q.r)
            edge_method = 'grid';
            [grid, density, m, v, mass_at_0, K_hat, l_hat, u_hat, x, f_hat] = ...
                compute_esd_ode(q.t, q.w, q.gamma,q.epsilon,edge_method,q.M);
        else
            [grid, density, m, v, mass_at_0, K_hat, l_hat, u_hat, x, f_hat] = ...
                compute_esd_ode_non_atomic(q.t, q.w, q.gamma, q.r,q.w_int,q.epsilon,q.M);
        end
        
    case {'fp'}
        if isempty(r)
            grid = compute_esd_ode(t, w, gamma,epsilon,M);
            [density,m,v] = compute_esd_fp(q.t, q.w, q.gamma,q.epsilon);
        else
            grid = compute_esd_ode_non_atomic(t, w, gamma, r,w_int,epsilon,M);
            [density,m,v] = compute_esd_fp_non_atomic(q.t,q.w,q.r,q.w_int,q.gamma,q.epsilon, grid);
        end
        
    case {'newton'}
        if isempty(r)
            newton_update = 'null';
            [grid, density, m, v, mass_at_0, K_hat, l_hat, u_hat, x, f_hat] =...
                compute_esd_newton(q.t, q.w, q.gamma, newton_update, q.epsilon,q.M);
        else
            error('Newton method currently not available with non-atomic measures\n','Newton_non_atomic');
        end
end

end
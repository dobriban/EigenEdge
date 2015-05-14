function pa = parseArgumentsSpectrode(varargin)
%parses arguments in varargin; specifically designed for Atomic
%validation is a binary - specifies if validation is attempted

%Must create parameters to calls of the form
% [grid, density, m, v, mass_at_0, K_hat, l_hat, u_hat, x, f_hat] = ...
%    compute_esd_ode_non_atomic(t, w, gamma, r,w_int,epsilon,M)

pa = inputParser;


%algorithms and parameters
defaultAlg = 'spectrode';
expectedAlgs = {'spectrode','fp','Newton'};
w = [];
r = [];
w_int = [];
epsilon = 1e-4;
M = floor(sqrt(1/epsilon))+3;

%set defaults
%required
addRequired(pa,'t');
addRequired(pa,'gamma');

%optional
addOptional(pa,'w',w,@isnumeric);
addOptional(pa,'r',r,@isnumeric);
addOptional(pa,'w_int',w_int,@isnumeric);
addOptional(pa,'epsilon',epsilon,@isnumeric);
addOptional(pa,'M',M,@isnumeric);

addParameter(pa,'alg',defaultAlg,@(x) any(validatestring(x,expectedAlgs)));
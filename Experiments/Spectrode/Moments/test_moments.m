%Test moment computations for limit ESD
%% Test correctness by computing theoretical and numerical solutions
% Null case where all eigenvalues are equal to 1
%first check that the density integrates to 1
p = 200;
gamma = 1/2;
epsi_array = [4,6,8];
t = 1;
for i=1:3;
    gamma_plus = (1+sqrt(gamma))^2;
    gamma_minus = (1-sqrt(gamma))^2;
    MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;
    epsilon = 10.^(-epsi_array(i));
    w = 1;
    fun = @(x) 1;
    integral_of_density(i) = esd_moment(t,w,gamma,fun,epsilon);
end

%% Null case where all eigenvalues are equal to 1
%mean
p = 200;
gamma = 1/2;
epsi_array = [4,6,8];
t = 1;
for i=1:3;
    gamma_plus = (1+sqrt(gamma))^2;
    gamma_minus = (1-sqrt(gamma))^2;
    MP_density = @(x) 1/(2*pi*gamma)* sqrt(max((gamma_plus-x).*(x-gamma_minus),0))./x;
    epsilon = 10.^(-epsi_array(i));
    w = 1;
    fun = @(x) x.^2;
    mean(i) = esd_moment(t,w,gamma,fun,epsilon);
end

%For the null MP law, the mean=1, variance = gamma
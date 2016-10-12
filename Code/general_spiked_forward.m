function [lambda,cos_right,cos_left,m,v] = general_spiked_forward(ell, t, w, gamma)
%Compute the characteristics of the general spiked model
%Using the spectrode method

%Input
% ell - population spike location (null: ell=0)
% t - population eigenvalues >=0
% w - mixture weights >=0;
% gamma - aspect ratio p/n

%Output
%lambda - asymptotic sample spike location
%cos_right,cos_left- squares of the right/left singular value angles
%m,v - Stieltjes transform and companion ST

ep = 1e-6;
[~, ~, ~, ~, ~, ~, ~, u_hat,~, ~,grid,v] =   compute_esd_ode(t, w, gamma,ep);
v = real(v);
m = 1/gamma*v- (1-1/gamma)./grid;
%b_sq_ind = min(find(grid>=max(u_hat))); %index of the upper edge
b_sq_ind = max(find(grid<=max(u_hat))); %index of the upper edge
b_sq = grid(b_sq_ind); %the smallest element above the the upper edge on the grid
D = grid.*m.*v;
PT=D(b_sq_ind);
%plot(grid,real(D))
%not increasing on last interval?
q = length(D);
D = D(b_sq_ind:q);
m = m(b_sq_ind:q);
v = v(b_sq_ind:q);
grid = grid(b_sq_ind:q);

v_p = @(v) (1/v^2 - gamma* sum( w.*(t.^(-1) + v).^(-2)))^(-1);
v_prime =  zeros(length(v),1);
for i=1:length(v)
    v_prime(i) = v_p(v(i));
end

%m = 1/gamma*v- (1-1/gamma)./grid;
m_prime =  1/gamma*v_prime+ (1-1/gamma)./(grid.^2);

%D = x*m*v
%D' = mv + x*m*v'+x*m'*v
D_prime = m.*v+grid.*(m_prime.*v+v_prime.*m);

if ell>1/PT
    lambda_ind = min(find(D<1/ell));
    lambda = grid(lambda_ind); %lambda = D^{-1}(1/ell)
    cos_right = m(lambda_ind)/(D_prime(lambda_ind)*ell); %m(lambda)/[D'(lambda)*ell];
    cos_left = v(lambda_ind)/(D_prime(lambda_ind)*ell);
else
    lambda = b_sq;
    cos_right = 0;
    cos_left = 0;
end

%vq = interp1(x,v,xq) interpolated values of a 1-D function using linear interpolation
v = interp1(grid,v,ell);
m = 1/gamma*v- (1-1/gamma)./ell;
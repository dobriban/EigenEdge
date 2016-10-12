function [lambda,cos_right,cos_left,m,v] = standard_spiked_forward(ell,gamma)
%Compute the characteristics of the standard spiked model
%Using explicit formulas

%Input
% ell - population spike location (null: ell=0, so that eigenvalue of population
% covariance matrix is 1+ell)
% gamma - aspect ratio p/n

%Output
%lambda - asymptotic sample spike location (eigenvalue)
%cos_right,cos_left - squares of the right/left singular value angles
%m,v - Stieltjes transform and companion ST

k  = length(ell);
lambda = zeros(k,1);
cos_right  = zeros(k,1);
cos_left  = zeros(k,1);
v  = zeros(k,1);
gamma_minus = (1-sqrt(gamma))^2;
gamma_plus = (1+sqrt(gamma))^2;

for i=1:k
    if (ell(i)<gamma^(1/2))&&(ell(i)>-gamma^(1/2)) %BBP PT
        lambda(i) = (1+gamma^(1/2))^2;
        cos_right(i) = 0;
        cos_left(i) = 0;
    else
        lambda(i) = (1+ell(i))*(1+gamma/ell(i));
        cos_right(i) = (1-gamma/ell(i)^2)/(1+gamma/ell(i));
        cos_left(i) = (1-gamma/ell(i)^2)/(1+1/ell(i));
    end
    x = lambda(i);
        im_mult = 1;
        if (x>gamma_minus)&&(x<gamma_plus)
            im_mult = 1i;
        end
        v(i) = 1/(2*x)*(-(1+x-gamma)+im_mult*(abs((1+x-gamma)^2-4*x))^(1/2));
end
m = 1/gamma*v- (1-1/gamma)./lambda;
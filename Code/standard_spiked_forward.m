function [lambda,cos_right,cos_left] = standard_spiked_forward(ell,gamma)
%Compute the characteristics of the standard spiked model
%Using explicit formulas

%Input
% ell - population spike location (null: ell=0, so that eigenvalue of population
% covariance matrix is 1+ell)
% gamma - aspect ratio p/n

%Output
%lambda - asymptotic sample spike location (eigenvalue)
%cos_right,cos_left - squares of the right/left singular value angles

k  = length(ell);
lambda = zeros(k,1);
cos_right  = zeros(k,1);
cos_left  = zeros(k,1);
%BBP PT
for i=1:k
    if ell(i)>gamma^(1/2)
        lambda(i) = (1+ell(i))*(1+gamma/ell(i));
        cos_right(i) = (1-gamma/ell(i)^2)/(1+gamma/ell(i));
        cos_left(i) = (1-gamma/ell(i)^2)/(1+1/ell(i));
    else
        lambda(i) = (1+gamma^(1/2))^2;
        cos_right(i) = 0;
        cos_left(i) = 0;
    end
end
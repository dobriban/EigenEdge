function [lambda,cos_right,cos_left] = standard_spiked_forward(gamma,ell)
%Compute the characteristics of the standard spiked model
%Using explicit formulas

%Input
% gamma - aspect ratio p/n
% ell - population spike location (null: ell=0, so that eigenvalue of population
% covariance matrix is 1+ell)

%Output
%lambda - asymptotic sample spike location (eigenvalue)
%cos_right,cos_left - squares of the right/left singular value angles

%BBP PT
if ell>gamma^(1/2) 
    lambda = (1+ell)*(1+gamma/ell);
    cos_right = (1-gamma/ell^2)/(1+gamma/ell);
    cos_left = (1-gamma/ell^2)/(1+1/ell);
else
    lambda = (1+gamma^(1/2))^2;
    cos_right = 0;
    cos_left = 0;
end

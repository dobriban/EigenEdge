function [ell,cos_right,cos_left] = standard_spiked_inverse(lambda,gamma)
%compute characteristics of the spiked model based on empirical spike

%Inputs
%lambda - empirical spike
%gamma - aspect ratio

%Outputs
%ell - de-biased shrinkage of spike; estimate of population spike
%  null location: ell=0
%cos_right,cos_left - squares of the right/left singular value angles

%formula for ell: from Donoho-Gavish-Johnstone paper
%Optimal Shrinkage of Eigenvalues
%in the Spiked Covariance Model
%http://arxiv.org/pdf/1311.0851v2.pdf
%formula 6.5 p 17

p = length(lambda);
ell  = zeros(p,1);
cos_right = zeros(p,1);
cos_left = zeros(p,1);
for i=1:p
    if lambda(i)>(1+sqrt(gamma))^2
        a = lambda(i)+1-gamma;
        ell(i)=1/2*(a+sqrt(a^2-4*lambda(i)))-1;
        [~,cos_right(i),cos_left(i)] = standard_spiked_forward(ell(i),gamma); 
    else
        ell(i)=0;
        cos_right(i) = 0;
        cos_left(i) = 0;
    end
end
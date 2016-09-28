function [ell,cos_right,cos_left] = general_spiked_estimate(lambda,eigs,gamma)
%Estimate the characteristics of the general spiked model
%Input
% lambda - empirical spike
% eigs - empirical eigenvalues >=0
% gamma - aspect ratio p/n

%Output
% ell - estimated population spike location
%cos_right,cos_left- estimates of squares of the right/left singular value angles

%first order estimates
m_hat = mean(1./(eigs-lambda)); %Stieltjes transform
%m = 1/gamma*v- (1-1/gamma)./grid;
v_hat = gamma*m_hat-(1-gamma)/lambda; %companion ST
D_hat = lambda*m_hat*v_hat; %D-transform
ell = 1/D_hat; %pop spike

%Derivatives
m_prime_hat = mean(1./(eigs-lambda).^2); %Stieltjes transform
v_prime_hat = gamma*m_prime_hat+(1-gamma)/lambda^2; %companion ST

%D = x*m*v
%D' = mv + x*m*v'+x*m'*v
D_prime_hat = m_hat*v_hat+lambda*(m_prime_hat*v_hat+v_prime_hat*m_hat);
cos_right = m_hat/(D_prime_hat*ell);
cos_left = v_hat/(D_prime_hat*ell);



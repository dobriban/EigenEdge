function [m,v] = SparseStieltjes2(zee,gamma,eps,tau)
%compute the Stieltjes transform of the fractional rank model
%with population distribution (1-eps)*delta_1 + eps*delta_tau
%Other Inputs
%zee - grid where Stietltjes transform should be computed
%gamma - p/n

%Outputs
%m - Stieltjes Transform at zee
%v- Gram Stieltjes Transform at zee

%in the solution of the equation, we assume gamma = p/n
%because we solve the equation in the form
%-1/v = z - gamma*[ (1-eps)/(1+v) + eps/(tau^{-1}+v)]
d3 = zee; %coeff of cubic term
d2 = zee.*(1./tau+1) - gamma + 1;
d1 = zee./tau - gamma.*((1-eps)./tau+eps) + 1./tau+1;
d0 = ones(1,length(d3))./tau ;
 
% get the unique solution with positive imaginary part
v = zeros(length(d0),1);
mp_roots = zeros(3,length(d0));
for i= 1:length(d0),
    coef = [d3(i) d2(i) d1(i) d0(i)];
    
    r = roots(coef);
    mp_roots(:,i) = r;
    if (any(imag(r) > 0)),
        mximr = max(imag(r));
        v(i) = r(imag(r)==mximr);
    else
        if (max(imag(r))==0)&&(min(imag(r))==0)
        v(i) = 0; 
        else
            v(i) = NaN;
        end
    end
end

m = 1/gamma*v+ (1-1/gamma)./zee;
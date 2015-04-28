function s = random_discrete_draw(t,w,p)
% p random draws from the discrete distribution
% H = sum_k w_k Delta_t_k
%
% Inputs
% t,w - parameters of atomic distribution: H = sum w_i*delta_{t_i}
% p - number of samples

s = zeros(p,1);
c = [0;cumsum(w)];
for i=1:p
    x = unifrnd(0,1);
    j = max(find(c<x));
    s(i) = t(j);
end
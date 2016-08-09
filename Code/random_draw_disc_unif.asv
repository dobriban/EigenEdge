function s = random_draw_disc_unif(t,w,r,w_int,p)
% p random draws from the discrete +uniform mixture distribution
% H = sum w(i)*delta_t(i) + sum w_int(i)*unif[r(i,1), r(i,2)]
% where
% delta_a  - point mass at a
% unif(a,b) - uniform distribution on [a,b]
% Used in several Spectrode examples.

%
% Inputs
% t,w - parameters of atomic distribution: H = sum w_i*delta_{t_i}
% p - number of samples
% r - intervals where the continuous part of the spectrum is supported,
%       real matrix of size Jx2, i-th row encodes the endpoints [a,b] of
%       the i-th continuous component of the spectrum
% w_int - weights of the continuous part of the spectrum


s = zeros(p,1);

if ~isempty(w_int)
    w_new = [w; w_int];
    c = [0;cumsum(w_new)];
else
    c = [0;cumsum(reshape(w,[length(w),1]))];
end
for i=1:p
    x = unifrnd(0,1);
    j = find(c<x,1,'last');
    % draw a discrete or a continuous random variable
    if j<=length(w)
        s(i) = t(j);
    else
        s(i) = unifrnd(r(j-length(w),1),r(j-length(w),2));
    end
end

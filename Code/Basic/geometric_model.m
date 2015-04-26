function [t,w]  = geometric_model(a,b,K)

t = b.^((0:K-1)');
w = a.^((0:K-1)');
w = w/sum(w);

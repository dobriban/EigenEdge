function [t,w]  = geometric_model(a,b,K)
%used in several Spectrode tests

t = b.^((0:K-1)');
w = a.^((0:K-1)');
w = w/sum(w);

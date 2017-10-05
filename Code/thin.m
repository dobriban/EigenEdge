function [t,w] = thin(t,max_spec)

w = ones(max_spec,1)/max_spec;
t = quantile(t,max_spec)';
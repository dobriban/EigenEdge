function [t, w]  = comb_model(K, t_min,t_max,gap)
%compute a population spectrum from the comb model
%Used in Spectrode examples

if ( K*(K-1)*gap /2 >=1)
  error('The gap needs to be small enough that K*(K-1)*gap /2 <1');
end
t = linspace(t_min,t_max,K)';
unnorm_weights = (1:1:K)'*gap;
w = unnorm_weights - sum(unnorm_weights)/K +1/K;

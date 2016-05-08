function [K,l,u,final_l,final_u] = support_identification(grid,density,thresh)
%support identification using thresholding
%Used in the Support identification experiment in the Spectrode paper. 
support_int = find(1/pi*imag(density)>thresh);
non_support_int = find(1/pi*imag(density)<=thresh);
K = 0;
u_index = 0;
l_index = [];
while u_index(K+1)<=length(grid)-1
  K = K+1;
  if isempty(min(support_int(support_int>u_index(K))))
    K = K-1;
    break
  else
    l_index = [l_index; min(support_int(support_int>u_index(K)))]; %length K
    u_index = [u_index; min(non_support_int(non_support_int>l_index(K)))-1];
  end
end

%retain only intervals of length >1
u_index = u_index(2:K+1);
final_l = [];
final_u = [];
for i=1:K
  if (l_index(i) < u_index(i)- 1)
    final_l = [final_l; l_index(i)];
    final_u = [final_u; u_index(i)];
  end
end

l = grid(final_l);
u = grid(final_u);
K = length(u);
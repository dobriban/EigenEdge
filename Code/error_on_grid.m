function [l_err,u_err] = error_on_grid(f,K_hat, l_hat, u_hat, K,ind_l,ind_u)
%compute bound on error on the grid f
%Used in the Support identification experiment in the Spectrode paper. 

K_err = abs(K-K_hat);

if (K_err>0)
  l_err = inf;
  u_err = inf;
else
  l_err_arr = zeros(K,1);
  u_err_arr = zeros(K,1);
  
  for i=1:K
    %find the indices of the hat estimators
    ind_l_hat = max(find(f==l_hat(i)));
    ind_u_hat = min(find(f==u_hat(i)));
    
    %find the smaller & larger indices -ie the width of the maximum
    %approximation error consistent with these approximations
    mini = min(ind_l_hat,ind_l(i));
    mini = max(mini-1,1);
    maxi = max(ind_l_hat,ind_l(i));
    maxi = min(length(f),maxi+1);
    l_err_arr(i) = abs(f(maxi)-f(mini));
    
    mini = min(ind_u_hat,ind_u(i));
    mini = max(mini-1,1);
    maxi = max(ind_u_hat,ind_u(i));
    maxi = min(length(f),maxi+1);
    u_err_arr(i) = abs(f(maxi)-f(mini));  
  end
  l_err = mean(l_err_arr);
  u_err = mean(u_err_arr);
end
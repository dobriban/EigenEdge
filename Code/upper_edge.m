function u= upper_edge(t, gamma)

p = length(t); 
max_spec = 100;
if p>max_spec
    t = thin(t,max_spec);
    p = max_spec;
end
v_L =-1/max(t);

%if epsilon is so large that the grid is empty, make it smaller
epsilon = 1e-4;
while (v_L+epsilon)>=(-epsilon)
    epsilon = epsilon/2;
end

w = ones(p,1)/p;
z = @(v) - 1/v + gamma* sum(w.*t./(1 + t.*v));
z_prime = @(v) 1/v^2 - gamma* sum(w.*t.^2./((1 + t.*v).^2));
if (z_prime(v_L+epsilon)>-Inf)&&(z_prime(v_L+epsilon)<Inf)
    n_it=0;
    while (z_prime(v_L+epsilon)*z_prime(-epsilon)>0)
        epsilon = epsilon/2;
        n_it = n_it+1;
    end
    v_s = fzero(z_prime,[v_L+epsilon,-epsilon]);
    u = z(v_s);
else
    u = Inf;
end


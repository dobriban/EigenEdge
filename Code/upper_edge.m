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
M = floor(sqrt(1/epsilon))+3;
v = linspace(v_L+epsilon,-epsilon, M); %adaptive grid-forming

w = ones(p,1)/p;
[~,u] = evaluate_inverse_ST(t,w,gamma, v);


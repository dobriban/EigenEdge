function y = esd_all_quantiles_grid(grid,density,mass_at_0,q)
%Compute a quantile of the limit spectrum of covariance matrices
%based on the pre-computed grid and density

%here q can be a vector of quantiles
grid = [0; grid];
l = length(grid);
pr = zeros(l,1);
pr(2:l) = (grid(2:l)-grid(1:l-1)).*density;
pr(1) = mass_at_0;
y = zeros(length(q),1);
y_sort = zeros(length(q),1);
ind=1;

cdf = cumsum(pr);
%q_sort(i)  = q(ind(i))
[q_sort, sort_ind] = sort(q);
%we assume q_sort is "smaller" than the length of the grid
%so it makes sense to iterate over cdf/grid inside the loop
for i=1:length(q_sort)
    while (ind<length(cdf))&&(cdf(ind)<q_sort(i))
        ind=ind+1;
    end
    y_sort(i) = grid(ind);
end 
%unsort 
for i=1:length(q)
    y(sort_ind(i)) = y_sort(i);
end

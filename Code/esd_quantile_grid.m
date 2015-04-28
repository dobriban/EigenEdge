function y = esd_quantile_grid(grid,density,mass_at_0,q)
%Compute a quantile of the limit spectrum of covariance matrices
%based on the pre-computed grid and density

grid = [0; grid];
l = length(grid);
pr = zeros(l,1);
pr(2:l) = (grid(2:l)-grid(1:l-1)).*density;
pr(1) = mass_at_0;

cdf = cumsum(pr);
%plot(grid,cdf);
ind_min=find(cdf>=q,1,'first');
y = grid(ind_min);
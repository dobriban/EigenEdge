function y = esd_moment_grid(grid,density,mass_at_0,fun)
%Compute a moment of the limit spectrum of covariance matrices
%based on the pre-computed grid and density
Y = fun(grid).*density;
y = trapz(grid,Y);

if (mass_at_0>0)
    y = mass_at_0*fun(0) + y;
end

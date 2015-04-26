function y = esd_mode_grid(grid,density)
%Compute the mode of the limit spectrum of covariance matrices
%based on the pre-computed grid and density

mode_ind = find(density==max(density));
y = grid(mode_ind);

function [n_points, z_grid, v_grid] =  grid_output(v_local,z)
%Used by Spectrode in compute_esd_ode, with Brent method

n_points = length(v_local); %num grid points where increasing
c = length(v_local);
for i=1:c
    z_grid(i) = z(v_local(i)); %z-space
end
v_grid = v_local; %in the v-space
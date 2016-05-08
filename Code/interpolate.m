function v_x = interpolate(x,grid,v_grid)
%interpolate the function "v_grid" defined on the grid "grid"
%onto the new grid "x"
%Used in computing optimal LSS

%Uses linear interpolation
j=1;
v_x = zeros(length(x),1);
for i=1:length(x)
    while (grid(j)<x(i))&&(j<length(grid))
        j = j+1;
    end %end up with j st grid(j)>=x(i)>grid(j-1)
    if j==1
        v_x(i) = v_grid(j);
    else %grid(j-1)<x(i)<=grid(j)
        epsi = (x(i) - grid(j-1))/ (grid(j) - grid(j-1));
        v_x(i) = epsi*v_grid(j-1)+(1-epsi)*v_grid(j);
    end
end
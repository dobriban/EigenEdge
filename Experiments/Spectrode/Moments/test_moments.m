%test the computation of moments of the standard MP law
%compare against true values; print to file
gamma = 1/2;
w = 1;
t = 1;
epsilon = 1e-8;
[grid, density,~, ~, mass_at_0] =  compute_esd_ode(t, w, gamma,epsilon);

%%
functions = {'x','\\log(x)','\\log^2(x)'};
%%
A = 3;
num_val = zeros(A,1);
true_val = zeros(A,1);
%% x
num_val(1) = esd_moment_grid(grid,density,mass_at_0,@(x)x);
true_val(1) = 1;

%% log x
num_val(2) = esd_moment_grid(grid,density,mass_at_0,@(x)log(x));
true_val(2) =  - 1 + (gamma - 1)/gamma*log(1-gamma);

%% x* log x
num_val(3) = esd_moment_grid(grid,density,mass_at_0,@(x)log(x).^2);
true_val(3) = NaN;

%%
rel_err = abs(true_val-num_val);

%% Latex
%print high level statistics
filename = ['test_moments_latex_output'];
fileID = fopen([filename '.txt'],'w');

fprintf(fileID,['Results of Computing Moments\n\n']);

for i=1:A
    str = sprintf('$%s$\t',functions{i});
    fprintf(fileID,str);
    fprintf(fileID,' & ');
    fprintf(fileID,num2str(true_val(i)));
    fprintf(fileID,' & ');
    fprintf(fileID,num2str(num_val(i)));
    fprintf(fileID,' & ');
    fprintf(fileID,num2str(rel_err(i)));
    fprintf(fileID,'\\\\ \\hline \n');
end


fclose(fileID);
fprintf(['Saved Results to ./Results/' filename '.txt\n']);

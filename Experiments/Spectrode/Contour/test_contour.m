%test the computation of contour integrals of the standard MP law
%compare against true values; print to file
gamma = 1/2;

%%
functions = {'x','\\log(x)','\\log^2(x)'};
%%
A = 3;
num_val = zeros(A,1);
true_val = zeros(A,1);
%%
g = @(x) x;
num_val(1) = LSS_mean_standard_mp(g, gamma);
true_val(1) = 0;
%%
g = @(x) log(x);
num_val(2) = LSS_mean_standard_mp(g, gamma);
true_val(2) = 1/2*log(1-gamma);
%%
g = @(x) log(x).^2;
num_val(3) = LSS_mean_standard_mp(g, gamma);
true_val(3) = NaN;

%%
err = abs(true_val-num_val);

%% Latex
%print high level statistics
filename = ['test_contour_latex_output'];
fileID = fopen([filename '.txt'],'w');

fprintf(fileID,['Results of Computing Contour Integrals\n\n']);

for i=1:A
    str = sprintf('$%s$\t',functions{i});
    fprintf(fileID,str);
    fprintf(fileID,' & ');
    fprintf(fileID,num2str(true_val(i)));
    fprintf(fileID,' & ');
    fprintf(fileID,num2str(num_val(i)));
    fprintf(fileID,' & ');
    fprintf(fileID,num2str(err(i)));
    fprintf(fileID,'\\\\ \\hline \n');
end


fclose(fileID);
fprintf(['Saved Results to ./Results/' filename '.txt\n']);
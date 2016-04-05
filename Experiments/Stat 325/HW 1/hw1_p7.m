%% Partial solution to Problem 7 in HW1, from STAT 325
cd('C:\Git\EigenEdge\Experiments\Stat 325\HW 1')
addpath('..\..\..\Code')
%% Compute quantiles of MP law
t = 1;
w = 1;
N = 100;
gamma_v = linspace(2,200,N);
p = 197146; %from slides
n = 1386; %1387, but lose 1 degree of freedom by centralization (from file)
gamma_null = p/n; %the forma value of gamma
prop_nonnull = 1 - 1/gamma_null;
prop_null = 1/gamma_null;
%%
%L  = 0.25; U  = 0.75;
L  = 0.35; U  = 0.65;
q1  = prop_nonnull + 0.25*prop_null;
q3  = prop_nonnull + 0.75*prop_null;
q1_mp = zeros(N,1);
q3_mp = zeros(N,1);
tic
for i=1:N
    str = sprintf('gamma %d out of %d.\n',i,N);
    fprintf(str);
    toc
    q1_mp(i) = esd_quantile(t, w, gamma_v(i), q1);
    q3_mp(i) = esd_quantile(t, w, gamma_v(i), q3);
end

%%
load('Eigs.mat','Eigs')
r = q3_mp./q1_mp;
r_data = quantile(Eigs,U)/quantile(Eigs,L);
%%
figure, hold on
plot(gamma_v,r);
plot(gamma_v,r_data*ones(N,1))
legend('R_G','location','Best')
xlabel('gamma')

%% save figures
filename = sprintf( 'R_G.png');
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% where does q1 become 0?
%1 - 1/g_pt = q1
g_pt = 1/(1-q1);
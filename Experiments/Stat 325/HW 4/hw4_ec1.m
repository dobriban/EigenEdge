%% Partial solution to Extra Credit Problem 1 of HW4, from STAT 325
%%Analytic sampling variance of the CDF of Wishart eigenvalues
cd('C:\Git\EigenEdge\Experiments\Stat 325\HW 4')
addpath('..\..\..\Code')
%% load eigenvalues
gamma  = 1/2;
load('hw4_p6_raw_data.mat')
%% Part 2.
%% Approximate the empirical variability in F_n(x), by that in Phi_n,h(x)
x = 1;
h_grid = 10.^(-8:1:0);
var_cdf = zeros(length(h_grid),1);
for i=1:length(h_grid)
    h=h_grid(i);
    phi = @(l) 1./(1 + exp( (l-x)./h));
    lss_phi = @(v) sum(phi(v));
    lss = zeros(n_m,1);
    for k=1:n_m
        lss(k) = lss_phi(eigs(k,:));
    end
    var_cdf_approx(i) = var(lss);
end

%% variance of F_n
cdf = sum(eigs' < x);
var_cdf_true = var(cdf);
%% plot variance and approximations
figure, hold on
h1 = plot(-log10(h_grid), var_cdf_approx,'LineWidth',2);
h2 = plot(-log10(h_grid), var_cdf_true*ones(length(h_grid),1),'LineWidth',2);

xlabel('-log_{10}(h)')
ylabel('Var[\Phi_{x,h}(q)]')
legend([h1 h2], {'Approx', 'True'}, 'location', 'best')
set(gca,'FontSize',14);
%%
filename = sprintf( './variance_of_MP_quantiles_approx_gamma=%.2f.png',gamma);
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% Part 3
%% Calculate the limit of the contour integral
x = 1;
h_grid = 10.^(-6:1:0);
r_grid = 1+10.^(-5:1:-1);
d = sqrt(gamma);
var_contour = zeros(length(h_grid),length(r_grid));
for k=1:length(r_grid)
    r = r_grid(k);
    for i=1:length(h_grid)
        h=h_grid(i);
        phi = @(l) 1./(1 + exp( (l-x)./h));
        
        ci = @(theta) cos(theta) + 1i*sin(theta);
        ci_p = @(theta) -sin(theta) + 1i*cos(theta);
        
        c = @(u,v) -1/(4*pi^2).*phi(abs(1+d.*ci(u)).^2).* phi(abs(1+d.*ci(v)).^2).*...
            ci_p(u).*ci_p(v)./(ci(u)-r*ci(v)).^2;
        %2D integral
        %q = integral2(fun,xmin,xmax,ymin,ymax)
        var_contour(i,k) =integral2(c,0,2*pi,0,2*pi);
    end
end
%% plot as a function of r and h
%diverges as r->1
figure, hold on
h1 = plot(-log10(h_grid(1:(length(h_grid))-1)), log10(real(var_contour(1:(length(h_grid)-1),:))),'LineWidth',2);

xlabel('-log_{10}h')
ylabel('log_{10}Var[\Phi_{x,h}(q)]')
legend({'r=1e-5', 'r=1e-4','r=1e-3','r=1e-2','r=1e-1'}, 'location', 'best')
set(gca,'FontSize',14);
%%
filename = sprintf( './variance_of_MP_quantiles_contour_int_approx_gamma=%.2f.png',gamma);
saveas(gcf, filename,'png');
fprintf(['Saved Results to ' filename '\n']);

%% Simpler example showing the ill-posedness of the contour integral
x = 1;
r= 1+1e-5;
d = sqrt(gamma);

ci = @(theta) cos(theta) + 1i*sin(theta);
ci_p = @(theta) -sin(theta) + 1i*cos(theta);

c = @(u,v) -1/(4*pi^2).*abs(1+d.*ci(u)).^2.*abs(1+d.*ci(v)).^2.*...
    ci_p(u).*ci_p(v)./(ci(u)-r*ci(v)).^2;
J_1 =integral2(c,0,2*pi,0,2*pi);
%Get 5*10^6
%% Conclusions
%Using J_1 formula numerically is a bad idea
%because of the division by (xi_1-r*xi_2)^2 - near singular
%the integral explodes for r near 1
%The limit itself exists, but this is a bad way to look at it

% it is better to use the equivalent formula in terms of the covariance
% kernel
%there, the limit r->1 has been evaluated already
%also, there is no more contour integral (there are some cancellations due
%to symmetry)

%% So, use the formula
% J_1(f,g) = \int f'(x) g'(y) k(x,y) dx dy

%1. Empirically, does J_1(phi_{x,h},phi_{x,h}) tend to a limit?
%2. Theoretically? 

K = cov_kernel_standard_mp(x,y,gamma);

function [grid_alt,LSS,density_alt,asy_effect_size,K_hat,l_hat,u_hat,upper_pt,ind_in,any_outside] = ...
    optimal_LSS(t,s_null,s_alt,gamma,LSS_comput_method,n,experimental_mode,epsi)
%compute the optimal LSS for testing against a spiked alternative
%H_0: spectrum = [t; s_null]; against
%H_1: spectrum = [t; s_alt];

%Inputs
%t - null eigenvalues
%s_null, s_alt - loc of null and alt spikes; assumes s_null =/= s_alt
%gamma - aspect ratio
%LSS_comput_method - computational method for LSS
%n - sample size (needed to set LSS above the PT)
%experimental_mode  - enable experimental mode (default = 1)
%epsi - accuracy for Spectrode

%Outputs
%grid - grid where LSS is computed
%LSS - the optimal LSS on grid
%density - density of the derivative of the MP map
%asy effect_size - the effect size of the test
% K_hat, l_hat, u_hat - SPECTRODE outputs (support of ESD)
% upper_pt - upper phase transition threshold

%a default choice for n
if ~exist('n','var')
    n = floor(length(t)/gamma);
end

%enables some experimental features
if ~exist('experimental_mode','var')
    experimental_mode = 1; 
end

%accuracy for Spectrode
if ~exist('epsi','var')
    epsi = 5*1e-6;
end

%epsi = 1e-4;
w = 1/length(t)*ones(length(t),1);
w_alt = 1/length(s_alt)*ones(length(s_alt),1);
[grid_alt, density_alt,~,ind_in,K_hat,l_hat,u_hat,a,b,upper_pt,v_grid_alt] = MP_derivative(t,w,s_alt,w_alt,gamma,epsi);

%in experimental mode, recompute the LSS if the spike is right above the PT
if experimental_mode == 1
  if (length(s_alt)==1)&&(s_alt>upper_pt)&&(s_alt<0.75*(1+sqrt(gamma))*upper_pt)
      s_alt = 0.99*upper_pt;
      [grid_alt, density_alt,~,ind_in,K_hat,l_hat,u_hat,a,b,upper_pt,v_grid_alt] = MP_derivative(t,w,s_alt,w_alt,gamma,epsi);
  end
end

LSS = zeros(length(grid_alt),1);
%plot(grid,density)
 
%If there is any spike outside the spectrum, then the asymptotic power is
%unity
r = -1./s_alt;
any_outside = 0;
ind_outside = zeros(length(r),1);
for i=1:length(r)
    for j=1:length(a)
        %the i-th pop spike satisfies the criterion for being outside the
        %bulk; in the interval [a(j), b(j)], which is in the image of the complement of S under v
        if (a(j)<r(i))&&(r(i)<b(j))
            any_outside = 1;
            ind_outside(i) = 1;
        end
    end
end

z = @(v) - 1./v + gamma* mean(t./(1 + t.*v));
psi_prime = @(alpha) 1-gamma*mean(t.^2./(alpha-t).^2);
%kernel function to use in LSS if we are above the PT
%Kern = @(x) max(0,1-abs(x));
%%Kern = @(x) min(1,max(0,x+1)); %not symm
%%Kern = @(x) exp(-x^2/2); discontinuous at edges
%Kern = @(x) min(1,max(0,2-abs(2*x)));
%Kern = @(x) min(1,max(0,3-abs(3*x)));
Kern = @(x) max(0,1-x^2);
%figure, ezplot(Kern)
%set the LSS to be 1 in the right neighborhoods
if any_outside==1
    for i=1:length(r)
        %if the i-th pop spike is pushed outside the bulk
        if (ind_outside(i) ==1)
            pop_s = s_alt(i);
            sample_s = z(-1/pop_s); %use the alternative form of the spike forward map
            %ind_plus = sum(grid_alt<sample_s); %find the exact location in the grid
            %set a neighborhood of the LSS around sample_s to be 1
            %how large a neighborhood?
            %intuitively, the fluctuations of the sample spike are or order
            %p^{-1/2} above the PT. And I would like to capture them all.
            %so should go c*p^{-1/2} for some decently large c
            %the spike variance for pop spike pop_s is (see Yao et al. book
            %2015, Thm 11.11, pg 235. )
            %2*pop_s^2*psi'(pop_s) = 2/ v'(sample_s)
            vari = 2*pop_s^2*psi_prime(pop_s);
            %3 standard errors are needed (at least) to get close to
            %full power for large values of the spike
            num_sd = 3; 
            width = num_sd*sqrt(vari)*n^(-1/2); 
            ind_plus = find(abs(grid_alt-sample_s)<=width);
            % deal with a corner case
            if ind_plus==0
                ind_plus=1;
            end
            for j=1:length(ind_plus)
                x = (grid_alt(ind_plus(j))-sample_s)/width;
                LSS(ind_plus(j))=Kern(x);  
            end
            %set LSS to 1 in the outermost regions
            if sample_s>max(u_hat)
                LSS(grid_alt>sample_s)=1;
            end
            if sample_s<min(l_hat)
                LSS(grid_alt<sample_s)=1;
            end
        end
    end
    asy_effect_size=Inf;
    
else
    %if the distributions are fully overlapping
    %compute the weak derivative wrt the null spikes
    w_null = 1/length(s_null)*ones(length(s_null),1);
    [grid_null, density_null] = MP_derivative(t,w,s_null,w_null,gamma,epsi);
    %interpolate onto the alternative grid
    density_null_int = interp1(grid_null, density_null,grid_alt);
    
    dgrid = grid_alt(2:length(grid_alt))-grid_alt(1:length(grid_alt)-1);
    dgrid = [dgrid; dgrid(length(dgrid))];
    density=density_alt-density_null_int;
    L = cumsum(density.*dgrid);
    
    %compute covariance kernel 
    if ~exist('LSS_comput_method','var')
        LSS_comput_method = 'collocation';
    end
    switch LSS_comput_method
        case 'collocation'
            %LSS_diff= compute_lss_colloc(grid_alt(ind_in),L(ind_in),t,w,gamma);
            LSS_diff= compute_lss_colloc_2(grid_alt(ind_in),L(ind_in),v_grid_alt(ind_in),t,w,gamma,epsi);
        case 'diag_regularization'
            K = cov_kernel(grid_alt(ind_in),grid_alt(ind_in),t,w,gamma);
            reg = 1e-4*trace(K)/sum(ind_in);
            K = K + reg*eye(sum(ind_in));
            
            %solve discretized equation
            LSS_diff = -K\L(ind_in);
    end
    LSS(ind_in) = cumsum(LSS_diff.*dgrid(ind_in));
    
    
    %grid_alt - main grid
    %L - distribution function of weak derivative (pointwise evaluation on grid)
    %LSS - optimal LSS (pointwise evaluation on grid)
    %LSS_diff - (K^+)*L
    asy_effect_size = -LSS_diff'*(L(ind_in).*dgrid(ind_in));
    asy_effect_size = length(s_alt)*asy_effect_size^(1/2);
    %sometimes it's imaginary
    
    %standardize
    LSS(ind_in) = LSS(ind_in)/max(abs(LSS));
    
    %extend by linear interpolation to the outside of the spectrum
    edges = ind_in(2:length(ind_in))-ind_in(1:length(ind_in)-1);
    l_edges = find(edges==1)+1; %lower edges of spectrum
    u_edges = find(edges==-1); %upper edges of spectrum
    l=length(u_edges);
    L = length(LSS);
    
    %extend LSS as constant on the two sides
    LSS(1:l_edges(1)) = LSS(l_edges(1));
    LSS(u_edges(l):L) = LSS(u_edges(l));
    
    %linear interpolation outside clusters
    for i=1:l-1
        x = [grid_alt(u_edges(i));grid_alt(l_edges(i+1))]; %sample points
        v = [LSS(u_edges(i));LSS(l_edges(i+1))]; %sample values
        xq = grid_alt(u_edges(i)+1:l_edges(i+1)-1); %query points
        LSS(u_edges(i)+1:l_edges(i+1)-1) = interp1(x,v,xq);
    end
    
end



function [X, TrueSpectrum] = generateCorrelatedData(n,p,testType,betaParams,diff)
%Generate population spectra and data

% Generate Data with specified population co
% Inputs
% n - sample size
% p - dimension
% other inputs -parameters of covariance matrix
% testType - the name of the spectrum, see list specified below
% betaParams - parameters of beta mixture
%
% Outputs
% X - n x p data matrix
% TrueSectrum - spectrum of population covariance matrix



switch testType
    
    case 'ones'
        TrueSpectrum = ones(p,1);
        
    case 'cluster'
        numClusters = 3;
        TrueSpectrum = floor(numClusters*rand(p,1))+1;
        
    case 'two-point'
        numClusters = 2;
        TrueSpectrum = ones(p,1)+diff*(floor(numClusters*rand(p,1)));
    case 'MP'
        
        X = 1/sqrt(n)*randn(n,p)*sqrt(Sigma);
        TrueSpectrum = svd(X,'econ').^2;
        
    case 'beta'
        TrueSpectrum = 1+10*betarnd(1,10,p,1);
        
    case 'toeplitz'
        a = 0.3;
        r = a.^(0:1:p-1);
        Sigma = toeplitz(r);
        TrueSpectrum = eig(Sigma);
        
    case 'bimodal'
        a = floor(p/2);
        b=p-a;
        TrueSpectrum=[1+betarnd(2,2,a,1); 3+ betarnd(2,2,b,1)];
        
    case 'beta_mixture'
        
        % 'beta_mixture'
        % Each component of the mixture is characterized by 5 parameters,
        % which are specified in the columns of the k by 5 betaParams array.
        % The 5 parameters (say in a vector \verb+pa+) have the following meaning:
        % the corresponding mixture is a shifted and scaled beta distribution:
        %
        %  pa(1)+pa(2)*beta(pa(3),pa(4))
        %
        %  with weight pa(5). Thus the sum of the 5th components of \verb+betaParams+
        %  must be 1. Then, the approximate number of elements from the given distribution
        %  is $\lfloor \verb+pa(5)+ \times p \rfloor$.
        %  Within each component of the mixture, the entries are iid draws
        %  from the specified beta distribution.
        
        TrueSpectrum = zeros(p,1);
        currentIndex = 0;
        betaNumMix = size(betaParams,1);
        for k=1:betaNumMix
            par= betaParams(k,:);
            L = floor(par(5)*p);
            TrueSpectrum(currentIndex+1:currentIndex+L) = par(1)+par(2)*betarnd(par(3),par(4),L,1);
            currentIndex = currentIndex+L;
        end
        if (currentIndex<p)  %fill in remainder
            r = p-currentIndex;
            TrueSpectrum(currentIndex+1:p) = par(1)+par(2)*betarnd(par(3),par(4),r,1);
        end
        
end
TrueSpectrum = sort(TrueSpectrum);

Sigma = diag(TrueSpectrum);
sqrtSigma = diag(sqrt(TrueSpectrum));
X = 1/sqrt(n)*randn(n,p)*sqrtSigma;

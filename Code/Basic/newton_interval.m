function [grid, v0] = newton_interval(t, w, gamma, endpoint_1, endpoint_2, M, m_start_point,m_end_point,newton_update)
%Compute limit spectrum of covariance matrices in a specific interval using Newton's method.
%Uses the equation in the form of Silverstein & Choi'95 (sometimes also
%called 'Silverstein equation')
%Method was obtained from Jack W. Silverstein  via personal communication.

%Written by Edgar Dobriban
% Inputs
% t - population eigenvalues >=0
% gamma - aspect ratio p/n
% w - mixture weights >=0; default: w = uniform
% endpoint_1, endpoint_2 - lower & upper endpoints where evaluate
% num - number of grid points
% start_point - real value of the starting point for the first newton
%       iteration (this will have been m(endpoint_1) )
% end_point - end point for the last newton
%       iteration (this will have been m(endpoint_2) )
% newton_update - method to find starting point for Newton method
%
%Outputs
%grid - real grid where Stieltjes transform is evaluated
% v0 - dual Stieltjes transform on real line

%Sanity check v(endpoint_1) = start_point
%equivalently v_inv(start_point) = endpoint_1
v_inv = @(v) - 1/v + gamma* sum(w.*t./(1 + t.*v));
err = abs(v_inv(m_start_point)-endpoint_1);
if (err>1e-4)
  fprintf('Mis-specified starting point in Newton Method\n');
end
err = abs(v_inv(m_end_point)-endpoint_2);
if (err>1e-4)
  fprintf('Mis-specified endpoint in Newton Method\n');
end

delta = 1e-3;
grid = linspace(endpoint_1 + delta, endpoint_2 -delta,M)';
%grid = linspace(endpoint_1, endpoint_2,M)';
v0 = zeros(M, 1);
spacing = grid(2)-grid(1);
max_restart = 0; %how many restarts allowed at maximum
%m starting point: add small imaginary component
%it appears that this parameter is crucial
imaginary_push = sqrt(delta);  %square root because it behaves like a square root function
start_point_i = m_start_point + 1i*imaginary_push;
epsi = 1e-4; %error tolerance
max_iter = 100;
rng(0);

for j=1:M
  x = grid(j);
  inv = @(v) - x - 1/v + gamma* sum(w.*t./(1 + t.*v));
  d_inv = @(v)  1/v^2 - gamma* sum(w.*t.^2./((1 + t.*v).^2));
  
  %first attempt
  for attempt  = 1:max_restart+1
    iter = 1;
    clear a;
    if attempt==1
      a(iter) = start_point_i;
    else %if diverged, try random restarts
      a(iter) = start_point_i + 1i*unifrnd(0,spacing);
    end
    while (abs(inv(a(iter)))>epsi) && (iter < max_iter)&& (abs(inv(a(iter)))<Inf)
      b = a(iter) - inv(a(iter))/d_inv(a(iter));
      a = [a b];
      iter = iter + 1;
    end
    
    if ~((abs(a(iter))==Inf)||(isnan(abs(a(iter))))||(imag(a(iter))==0))
      break
    end
  end
  %if converged to the conjugate solution, correct it
  if imag(a(iter))<0
    a(iter) = conj(a(iter));
  end
  v0(j) = a(iter);
  %this is the key update rule
  
  switch newton_update
    case 'null'
      start_point_i = a(iter);
    case 'lin_spacing'
      start_point_i = a(iter)+ 1i*spacing;
    case 'sqrt_spacing'
      start_point_i = a(iter)+ 1i*sqrt(spacing);
  end
end

%plot(f0,imag(m0))

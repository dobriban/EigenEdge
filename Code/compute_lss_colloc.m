function LSS = compute_lss_colloc_2(x,L,vx,t,w_null,gamma,epsi)
%compute the LSS by solving linear integral equation
%using collocation method
%New version

%Inputs
%gr - grid within support where to compute LSS
%L - vector containing pointwise evaluation of distribution function of
%   weak derivative on grid gr
%vx - values of v on grid
%t,w_null,gamma - the parameters of the null H
%epsi - accuracy in original grid

%Outputs
%L - LSS computed on downsampled gr

%compute v on a new denser grid
%epsilon = 1e-1*epsi; %this leads to inaccurate results for small epsi
epsilon = 1e-2*epsi;
[y, ~, ~, v] = compute_esd_ode(t, w_null, gamma,epsilon);

N=length(x);
%what index in the dense grid is closest to each index in the coarse grid?
ind = zeros(N,1);
for t=1:N
    [~, ind(t)] = min(abs(y - x(t)));
end

%vq = interp1(x,v,xq) returns interpolated values of a 1-D function at specific query points using linear interpolation.
%Vector x contains the sample points, and v contains the corresponding values, v(x).
%Vector xq contains the coordinates of the query points.

%interpolate v into x
%vx = interp1(y,v,x);

multi = 1.5;
W = zeros(N,N);
for t=1:N
    %vt = v(ind(t));
    vt = vx(t);
    for i=t:N
        int1 = 0;
        int2 = 0;
        %integrate on x(i-1), x(i)
        if i>1
            vi = v(ind(i-1):ind(i));
            K0 = 1/(2*pi^2)*log(1 + 4*imag(vi).*imag(vt)./abs(vi-vt).^2);
            if abs(i-t)<=1
                bad = isinf(K0);
                m = max(K0(~bad));
                K0(bad) = multi*m;
            end
            x0 = y(ind(i-1):ind(i))-x(i-1);
            dx = y(ind(i-1)+1:ind(i))- y(ind(i-1):ind(i)-1);
            dx = [dx(1); dx];
            int1 = sum(K0.*x0.*dx);
            int1 = 1/(x(i)-x(i-1))*int1;
        end
        %integrate on x(i), x(i+1)
        if i<N
            vi = v(ind(i):ind(i+1));
            K0 = 1/(2*pi^2)*log(1 + 4*imag(vi).*imag(vt)./abs(vi-vt).^2);
            if abs(i-t)<=1
                bad = isinf(K0);
                m = max(K0(~bad));
                K0(bad) = multi*m;
            end
            x0 = x(i+1)-y(ind(i):ind(i+1));
            dx = y(ind(i)+1:ind(i+1))- y(ind(i):ind(i+1)-1);
            dx = [dx(1); dx];
            int2 = sum(K0.*x0.*dx);
            int2 = 1/(x(i+1)-x(i))*int2;
        end
        W(t,i) = int1+int2;
    end
end

W = W + W' - diag(diag(W));
LSS= -W\L;
LSS(1)=0;


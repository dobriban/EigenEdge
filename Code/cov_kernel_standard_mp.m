function K = cov_kernel_standard_mp(x,y,gamma)
%compute the numerical values of the covariance kernel K controlling the
%variance of a LSS, in the standard MP model
%Note: heuristic multiplication on the diagonal

%Inputs
%x,y -real grids
%gamma - aspect ratio

%Outputs
%K - a matrix of dimension length(x) x length(y), where
%K(i,j) = k(x(i),x(j))
%and k(x,y) = 1/(2*pi^2)*log(1 + 4.*imag(v(x)).*imag(v(y))./abs(v(x)-v(y)).^2);
%v(x) is the companion ST of the MP distribution with aspect ratio gamma

K = zeros(length(x),length(y));
v_x = v(x,gamma);
v_y = v(y,gamma);
K_f(i,j)= 1/(2*pi^2)*log(1 + 4.*imag(v_x(i)).*imag(v_y(j))./abs(v_x(i)-v_y(j)).^2);
multi = 1.5;
for i = 1:length(x)
    for j = 1:length(y)
        if v_x(i)==v_x(j)
            if j>1
                K(i,j) = multi *K_f(i,j-1);
            else if i>1
                    K(i,j) = multi *K_f(i-1,j);
                end
            end
        else
            K(i,j) = K_f(i,j);
        end 
    end
end
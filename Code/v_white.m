function companion_st =v_white(z,gamma)
%compute the companion Stieltjes transform of white MP law
%evaluated on a grid z

%Input
%z - grid
%gamma - aspect ratio

%Output
%companion_st -  the companion Stieltjes transform of white MP law on z

gamma_minus = (1-sqrt(gamma))^2;
gamma_plus = (1+sqrt(gamma))^2;
companion_st = zeros(length(z),1);
for i=1:length(z)
    %two cases: real, or has nonzero imaginary part
    %real z(i): do it separately within an outside of the support
    if (imag(z(i))==0)
        x = z(i);
        im_mult = 1;
        if (x>gamma_minus)&&(x<gamma_plus)
            im_mult = 1i;
        end
        companion_st(i) = 1/(2*x)*(-(1+x-gamma)+im_mult*(abs((1+x-gamma)^2-4*x))^(1/2));
    else
        r = roots([z(i) z(i)+1-gamma 1]);
        p = imag(r);
        y = r(p>0);
        companion_st(i) = y(1);
    end
end
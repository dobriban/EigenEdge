function x =v(z,gamma)
%compute the companion Stieltjes transform of white MP law
x = zeros(length(z),1);
for i=1:length(z)
  r = roots([z(i) z(i)+1-gamma 1]);
  p = imag(r);
  y = r(p>0);
  x(i) = y(1);
end
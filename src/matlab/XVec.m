function x = XVec(d, j, n)

% XVec(d, j, n) returns a row vector of n values per interval for evaluating
% B-spline scaling functions of degree d, level j.

u = [0:n-1]/(n-1);
x = zeros(1,n*2^j);
p = 1;
for inter = 1:2^j,
  x(p:p+n-1) = (u + inter - 1)/2^j;
  p = p + n;
end;
return;

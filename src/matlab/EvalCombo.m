function f = EvalCombo(d, j, c, n)

% EvalCombo(d, j, c, n) evaluates a linear combination of B-spline scaling
% functions of degree d, level j, using coefficients c, at n points per
% interval.

w = BernsteinWeights(d, j) * c;
u = [0:n-1]/(n-1);
f = zeros(1, n*2^j);
p = 1;
for inter = 1:2^j
  for i = 0:d
    f(p:p+n-1) = f(p:p+n-1) + w(i+1+(d+1)*(inter-1)) * Bernstein(d, i, u);
  end;
  p = p + n;
end;
return;

function w = BernsteinWeights(d, j)

% w = BernsteinWeights(d, j) returns a matrix of B-spline scaling function
% weights for Bernstein polynomials of degree d, level j.

w = eye(2^j + d);
if d == 0
  return;
end;
u = Knots(d, j);
g = Greville(d, u);

for i = 1:2^j-1
  for r = 1:d
    [u, g, w] = InsertKnot(d, u, g, w, i/2^j);
  end;
end;
return;

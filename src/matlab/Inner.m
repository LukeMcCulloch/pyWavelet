function I = Inner(d, j)

% I = Inner(d, j) returns the inner product matrix for B-spline scaling
% functions of degree d at level j.

I0 = BernsteinInner(d);
n = 2^j + d;
I = zeros(n);
w = BernsteinWeights(d, j);
for k = 1:n
  w1 = reshape(w(:,k), d+1, 2^j);
  for l = 1:n
    w2 = reshape(w(:,l), d+1, 2^j);
    I(k,l) = trace(w1'*I0*w2);
    I(l,k) = I(k,l);
  end;
end;
I = I / 2^j;
return;

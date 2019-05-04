function f = Factorial(m)

% f = Factorial(m) returns the matrix of factorials of entries of m.

[r,c] = size(m);
f = zeros(r, c);
for i = 1:r
  for j = 1:c
    f(i,j) = prod(2:m(i,j));
  end;
end;
return;

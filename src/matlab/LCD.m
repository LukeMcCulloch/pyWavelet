function d = LCD(m)

% d = LCD(m) returns the least common denominator of the entries of m.

[num,denom] = rat(m, 1e-10);
[m, n] = size(denom);
d = 1;
for i = 1:m
  for j = 1:n
    g = gcd(d, denom(i,j));
    d = d*denom(i,j)/g;
  end;
end;
return;

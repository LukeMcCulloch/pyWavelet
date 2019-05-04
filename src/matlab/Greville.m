function x = Greville(d, u)

% x = Greville(d, u) returns the vector of Greville abscissa values
% corresponding to degree d and knot vector u.

l = length(u);
x = u(1:l-d+1);
for k = 2:d
  x = x + u(k:l-d+k);
end;
x = x / d;
return;

function P = FindP(d, j)

% P = FindP(d, j) returns the P matrix for B-spline scaling functions of
% degree d, level j.

d = fix(d);
if d < 0,
  error('FindP: Must have d >= 0.');
end;
j = fix(j);
if j < 1
  error('FindP: Must have j >= 1.');
end;

if d == 0
  P = [1; 1];
  for i = 2:j
    P = [P zeros(size(P)); zeros(size(P)) P];
  end;
else
  u = Knots(d, j - 1);
  g = Greville(d, u);
  P = eye(2^(j-1) + d);
  for k = 0:2^(j-1)-1
    [u, g, P] = InsertKnot(d, u, g, P, (2*k+1)/2^j);
  end;
end;
return;

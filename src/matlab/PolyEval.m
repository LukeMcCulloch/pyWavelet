function pret = PolyEval(g, p, gnew)

% pret = PolyEval(g, p, gnew) returns the values of a control polygon
% defined by abscissas g and ordinates p, evaluated at gnew.

[m, n] = size(p);
if length(g) ~= m
  error('PolyEval: Length of g and rows of p must be the same.');
end;

for i = 1:length(gnew)
  row = max(find(g <= gnew(i)));
  if row == m
    pret(i,:) = p(m,:);
  else
    frac = (g(row+1) - gnew(i))/(g(row+1) - g(row));
    pret(i,:) = frac*p(row,:) + (1 - frac)*p(row+1,:);
  end;
end;
return;

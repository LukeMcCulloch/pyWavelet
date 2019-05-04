function I = BernsteinInner(d)

% I = BernsteinInner(d) returns the matrix of inner products of Bernstein
% polynomials of degree d.

i = ones(d+1, 1)*[0:d];
j = i';
I = Choose(d, i).*Choose(d, j)./(Choose(2*d, i+j)*(2*d + 1));
return;

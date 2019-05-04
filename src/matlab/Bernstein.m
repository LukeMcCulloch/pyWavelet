function b = Bernstein(d, i, u)

% b = Bernstein(d, i, u) returns the value of the i'th Bernstein polynomial
% of degree d at u.  Here d >= 0, 0 <= i <= d, and 0 <= u <= 1.

b = Choose(d, i)*u.^i.*(1-u).^(d-i);
return;

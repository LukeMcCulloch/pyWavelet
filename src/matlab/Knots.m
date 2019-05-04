function x = Knots(d, j)

% x = Knots(d, j) returns a vector of knot values for B-spline scaling
% functions of degree d, level j.

x = [zeros(1, d-1) [0:2^j-1]/2^j ones(1,d)];
return;

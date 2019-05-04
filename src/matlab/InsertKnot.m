function [uret, gret, pret] = InsertKnot(d, u, g, p, unew)

% [uret, gret, pret] = InsertKnot(d, u, g, p, unew) inserts a new knot at
% unew for B-spline scaling functions of degree d, thereby modifying knot
% vector u, Greville abscissas g, and synthesis matrix p.

uret = sort([u unew]);
gret = Greville(d, uret);
pret = PolyEval(g, p, gret);
return;

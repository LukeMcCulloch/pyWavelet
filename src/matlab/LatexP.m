function P2 = LatexP(d, j, pattern)

% P2 = LatexP(d, j, pattern) prints out P for degree d, level j, pulling
% fraction out front.  Returns the actual matrix printed.  If pattern is
% present, LatexPattern is used to show the general pattern.

P = FindP(d, j);
denom = LCD(P);
P = round(P*denom);
denom = sprintf('%.0f', denom);
if ~strcmp(denom, '1')
  disp(['  \frac{1}{' denom '} ']);
end;
if nargin > 2
  LatexPattern(d, P, d + 2);
else
  P2 = Latex(P)/str2num(denom);
end;
return;

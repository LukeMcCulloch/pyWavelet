% MakeFigs creates a number of the figures that appear in "Wavelets for
% Computer Graphics: A Primer" by Eric J. Stollnitz, Tony D. DeRose, and
% David H. Salesin.

% One-dimensional Haar sequence of level-of-detail approximations.
HaarApprox('max', 'save');

% One-dimensional Haar sequence of compressed approximations.
HaarCompress('L2', 'save');

% B-spline scaling functions of degree d = 0 through 3, at level j = 1.
PlotScaling(3, 1, 'save');

% B-spline scaling functions and wavelets with d = 0 through 3, j = 3.
PlotAll(0, 3, 'save');
PlotAll(1, 3, 'save');
PlotAll(2, 3, 'save');
PlotAll(3, 3, 'save');
return;

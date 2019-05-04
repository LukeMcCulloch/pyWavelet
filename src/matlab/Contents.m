% Endpoint-interpolating B-spline wavelets.
% Eric J. Stollnitz (stoll@amath.washington.edu)
% Contents updated December 22, 1997.
% 
% Matlab functions available from:
%
%    http://www.cs.washington.edu/research/graphics/projects/wavelets/
%
% More information can be found in:
%
%    Wavelets for Computer Graphics: Theory and Applications.
%    Eric J. Stollnitz, Tony D. DeRose, and David H. Salesin.
%    Morgan Kaufmann, San Francisco, 1996.  ISBN 1-55860-375-1.
%
%    Wavelets for Computer Graphics: A Primer.
%    Eric J. Stollnitz, Tony D. DeRose, and David H. Salesin.
%    IEEE Computer Graphics and Applications, 15(3):76-84, May 1995 (part1),
%    and 15(4):75-85, July 1995 (part 2).
%
%
% Functions used to generate matrices and figures in the text
%   LatexMatrices
%   MakeFigs
% 
% Arbitrary degree endpoint-interpolating B-spline wavelets
%   FindP
%   FindQ
%   Inner
% 
% Haar basis transform and example functions
%   HaarApprox
%   HaarCompress
%   HaarCoords
%   HaarDecompStep
%   HaarDecomposition
%   HaarReconStep
%   HaarReconstruction
% 
% Input/output and error checking support functions
%   CheckInt
%   CheckLatex
%   Latex
%   LatexP
%   LatexPattern
%   LatexQ
%   PlotAll
%   PlotScaling
% 
% Arbitrary degree B-spline wavelet support functions
%   Bernstein
%   BernsteinInner
%   BernsteinWeights
%   Choose
%   EvalCombo
%   Factorial
%   Greville
%   InsertKnot
%   Knots
%   LCD
%   PolyEval
%   XVec

function [x, y] = HaarCoords(c, normalization)

% [x, y] = HaarCoords(c) returns x and y coordinates suitable for plotting a
% piecewise-constant function on [0,1] with scaling function coefficients c.
% Normalization can be either 'max' or 'L2'.

% Check length of c.
len = length(c);
j = log(len)/log(2);
if len ~= 2^j
  error('HaarCoords: c must be a vector of length 2^j.');
end;

% Scale appropriately for normalization.
if strcmp(normalization, 'L2')
  c = c*2^(j/2);
end;

% Construct x and y positions.
x = floor([0.5:0.5:len])/len;
y(1:2:2*len) = c;
y(2:2:2*len) = c;
return;

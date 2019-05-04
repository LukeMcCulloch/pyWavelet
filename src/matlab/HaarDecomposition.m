function y = HaarDecomposition(x, normalization)

% HaarDecomposition(x) returns the Haar wavelet decomposition of x.
% The normalization can be either 'max' or 'L2'.

% Check length of x.
len = length(x);
j = log(len)/log(2);
if len ~= 2^j
  error('HaarDecomposition: x must be a vector of length 2^j.');
end;

% Copy input, scaled appropriately for normalization.
if strcmp(normalization, 'L2')
  y = x/2^(j/2);
else
  y = x;
end;

% Repeatedly do summing and differencing on shorter and shorter part of y.
while len >= 2
  y = HaarDecompStep(y, len, normalization);
  len = len/2;
end;
return;

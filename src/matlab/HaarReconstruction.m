function y = HaarReconstruction(x, normalization)

% y = HaarReconstruction(x) returns the scaling function reconstruction of
% the Haar wavelet decomposition x.  The normalization can be either 'max'
% or 'L2'.

% Check length of x.
len = length(x);
j = log(len)/log(2);
if len ~= 2^j
  error('HaarReconstruction: x must be a vector of length 2^j.');
end;

% Repeatedly do summing and differencing on shorter and shorter part of y.
len = 2;
while len <= 2^j
  x = HaarReconStep(x, len, normalization);
  len = len*2;
end;

% Copy output, scaled appropriately for normalization.
if strcmp(normalization, 'L2')
  y = x*2^(j/2);
else
  y = x;
end;
return;

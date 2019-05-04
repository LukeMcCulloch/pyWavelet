function y = HaarReconStep(x, len, normalization)

% y = HaarReconStep(x, len, normalization) does one step of a Haar wavelet
% reconstruction on the first len elements of x, with normalization either
% 'max' or 'L2'.

if strcmp(normalization, 'L2')
  factor = sqrt(2)/2;
else
  factor = 1;
end;

y = x;
y(1:2:len-1) = (x(1:len/2) + x(len/2+1:len))*factor;
y(2:2:len) = (x(1:len/2) - x(len/2+1:len))*factor;
return;

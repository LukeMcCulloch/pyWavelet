function y = HaarDecompStep(x, len, normalization)

% y = HaarDecompStep(x, len, normalization) does one step of a Haar wavelet
% decomposition on the first len elements of x, with normalization either
% 'max' or 'L2'.

if strcmp(normalization, 'L2')
  factor = sqrt(2);
else
  factor = 2;
end;

y = [(x(1:2:len-1) + x(2:2:len))/factor ...     % sums
	(x(1:2:len-1) - x(2:2:len))/factor ...  % differences
	x(len+1:length(x))];                    % anything beyond len
return;

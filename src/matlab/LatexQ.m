function Q2 = LatexQ(d, j, pattern)

% Q2 = LatexQ(d, j, pattern) prints out Q for degree d, level j, showing
% normalization factor first if possible.  Returns the actual matrix
% printed.  If pattern is present, LatexPattern is used to show the general
% pattern.

% Find the Q matrix.
Q = FindQ(d, j, 'min');

% Find the inner product of each wavelet with itself.
I = Inner(d, j);
I_denom = LCD(I);
I = round(I*I_denom);
ip = diag(Q'*I*Q);

if d == 3 & j == 3
  % This case is annoying because none of the columns can be made integers.
  Q = Q*diag(1./sqrt(ip/I_denom));
  Q2 = Latex(Q);
else
  % Scale Q by normalization of middle columns.
  numer = I_denom;
  middle = ceil(length(ip)/2);
  denom = ip(middle);
  ip = ip/denom;
  Q = Q*diag(1./sqrt(ip));
  % Scale up by least common denominator to get integers.
  if any(round(Q(:,middle)) ~= Q(:,middle))
    lcd = LCD(Q(:,middle));
    Q = Q*lcd;
    denom = denom*lcd*lcd;
  end;
  % Make sure denominator really is an integer.
  if abs(round(denom) - denom) < 1e-4
    denom = round(denom);
  end;
  % If we're looking for a pattern, we need to take out 2^j.
  if nargin > 2
    numer = numer/2^j;
  end;
  % Simplify the fraction numer/denom.
  cd = gcd(numer, denom);
  numer = numer/cd;
  denom = denom/cd;
  % Convert to strings.
  numer = sprintf('%.0f', numer);
  denom = sprintf('%.0f', denom);
  % Print out normalization: sqrt(numer/denom).
  if ~strcmp(numer, '1') | ~strcmp(denom, '1') | nargin > 2
    s = ['  \sqrt{'];
    if ~strcmp(denom, '1')
      s = [s '\tinyfrac{'];
    end;
    if ~strcmp(numer, '1') | nargin <= 2
      s = [s numer];
    end;
    if ~strcmp(numer, '1') & nargin > 2
      s = [s ' \cdot '];
    end;
    if nargin > 2
      s = [s '2^j'];
    end;
    if ~strcmp(denom, '1')
      s = [s '}{' denom '}'];
    end;
    disp([s '}']);
  end;
  if nargin > 2
    % Print out pattern of Q.
    LatexPattern(d, Q, 3*d + 2);
  else
    % Print out matrix and get return value.
    Q2 = sqrt(str2num(numer)/str2num(denom))*Latex(Q);
  end;
end;
return;

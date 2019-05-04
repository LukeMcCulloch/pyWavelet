function Q = FindQ(d, j, normalization)

% Q = FindQ(d, j, normalization) 
%  returns the Q matrix for B-spline scaling functions 
%  of degree d, level j.  
%
%  If normalization is 'min' (or is not specified) 
%  then the smallest entry in each column is made 1.  
%
%  If normalization is 'max' 
%  then the largest entry in each column is made 1.
%
%  If normalization is 'L2' 
%  then the L_2 norm of each wavelet is made 1.

if nargin < 3
  normalization = 'min';
elseif ~strcmp(normalization, 'min') & ~strcmp(normalization, 'max') ...
      & ~strcmp(normalization, 'L2') & ~strcmp(normalization, 'lcd')
  error('FindQ: normalization must be ''min'', ''max'', or ''L2''.');
end;

P = FindP(d, j);
I = Inner(d, j);
if ~strcmp(normalization, 'L2')
%  P = round(P*LCD(P));
%  I = round(I*LCD(I));
end;

M = P'*I;
[m1, m2] = size(M);
n = m2 - rank(M);
Q = zeros(m2, n);
found = 0;
start_col = 0;
while (found < n/2) & (start_col < m2)
  start_col = start_col + 1 + (found > d);
  width = 0;
  rank_def = 0;
  while ~rank_def & (width < m2 - start_col + 1)
    width = width + 1;
    submatrix = M(:,start_col:start_col+width-1);
    rank_def = width - rank(submatrix);
  end;
  if rank_def
    % find null space of submatrix (should be just one column)
    q_col = null(submatrix);
    
    if strcmp(normalization, 'min')
      % normalize column so smallest nonzero entry has magnitude 1
      q_col = q_col/min(abs(q_col + 1e38*(abs(q_col) < 1e-10)));
    elseif strcmp(normalization, 'max')
      % normalize column so largest entry has magnitude 1
      q_col = q_col/max(abs(q_col));
    elseif strcmp(normalization, 'lcd')
      % normalize column to get all integers
      q_col = q_col/min(abs(q_col + 1e38*(abs(q_col) < 1e-10)));
      q_col = q_col*LCD(q_col);
    end;

    % change sign to give consistent orientation
    q_col = q_col*(-1)^(start_col + floor((d+1)/2) + (q_col(1,1) > 0));
    disp(q_col) %******TLM
    % correct any slight error for answers that should be integers
    if all(abs(submatrix*round(q_col)) < 1e-10) & any(round(q_col) ~= 0)
      q_col = round(q_col);
    end;
    
    % put column into strcmpleft half of Q
    found = found + 1;
    Q(start_col:start_col+width-1,found) = q_col;
    
    % use symmetry to put column into right half of Q in reverse order
    % and negated if degree is even
    Q(:,n-found+1) = flipud(Q(:,found))*(-1)^(d+1);
  end;
  disp (rank_def) 
  disp (width < m2 - start_col +1)
  disp (width)
  disp (m2)
  disp (start_col)
end;
  disp ('rank_def has kicked out') 
  disp (rank_def) 
if strcmp(normalization, 'L2')
  % normalize matrix so each column has L_2 norm of 1
  ip = Q'*I*Q;
  Q = Q*diag(1./sqrt(diag(ip)));
end;
return;

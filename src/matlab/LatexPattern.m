function LatexPattern(d, A, support)

% LatexPattern(d, A, support) outputs matrix A in LaTeX form, showing a
% pattern for degree d B-splines, with support given.

rows = (d + ceil(support/2))*2 + 5;
cols = 2*d + 6;
[m, n] = size(A);
disp('  \left[');
disp(['  \begin{array}{*{' num2str(d+3) '}{r}\at{\,}r\at{\,}*{' ...
	num2str(d+2) '}{r}}']);

% Figure out whether or not column contains only integers.
int = ones(1, n);
for j = 1:n
  for i = 1:m
    if abs(A(i,j) - round(A(i,j))) >= 5e-7
      int(j) = 0;
    end;
  end;
end;

% Figure out how wide each column should be.
width = zeros(1, n);
for i = 1:m
  for j = 1:n
    if A(i,j) ~= 0
      if int(j)
	len = length(sprintf('%g', A(i,j)));
      else
	len = length(sprintf('%.6f', A(i,j)));
      end;
      if len > width(j)
	width(j) = len;
      end;
    end;
  end;
end;

% Print out matrix.
for i = 1:rows
  s = '    ';
  for j = 1:cols
    if j <= d + 2
      if A(i,j) == 0
	s2 = blanks(width(j));
      elseif int(j)
	s2 = sprintf(['%' num2str(width(j)) 'g'], A(i,j));
      else
	s2 = sprintf(['%' num2str(width(j)) '.6f'], A(i,j));
      end;
    elseif j <= d + 5
      if i - d - ceil(support/2) == j - d
	s2 = '\cdot';
      else
	s2 = '     ';
      end;
    else
      ii = i - rows + m;
      jj = j - cols + n;
      if A(ii,jj) == 0
	s2 = blanks(width(jj));
      elseif int(jj)
	s2 = sprintf(['%' num2str(width(jj)) 'g'], A(ii,jj));
      else
	s2 = sprintf(['%' num2str(width(jj)) '.6f'], A(ii,jj));
      end;
    end;
    s = [s s2];
    if j < cols
      s = [s ' & '];
    else
      s = [s ' \\'];
    end;
  end;
  disp(s);
end;
disp('  \end{array} \right]');

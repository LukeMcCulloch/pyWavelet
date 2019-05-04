function A2 = Latex(A)

% A2 = Latex(A) outputs matrix A in LaTeX form, leaving out zeros,
% and returns the actual matrix that was output.

[m, n] = size(A);
disp('  \left[');
disp(['  \begin{array}{*{' num2str(n) '}{r}}']);

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

% Print out matrix and compute return values.
A2 = zeros(size(A));
for i = 1:m
  s = '    ';
  for j = 1:n
    if A(i,j) == 0
     s = [s blanks(width(j))];
    else
      if int(j)
	s2 = sprintf(['%' num2str(width(j)) 'g'], A(i,j));
      else
	s2 = sprintf(['%' num2str(width(j)) '.6f'], A(i,j));
      end;
      A2(i,j) = str2num(s2);
      s = [s s2];
    end;
    if j < n
      s = [s ' & '];
    else
      s = [s ' \\'];
    end;
  end;
  disp(s);
end;
disp('  \end{array} \right]');
return;

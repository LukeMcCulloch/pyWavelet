function PlotScaling(dmax, j, save)

% PlotScaling(dmax, j, save) plots B-spline scaling functions of degrees 0
% through dmax at level j.  If save is present, the plot is saved in an EPS
% file.

figure(1);
clg;
hold on;
j = j + 1;
n = 20;         % samples per interval
hspacing = 1.4; % horizontal spacing between functions
hsize = 0.8;    % horizontal size of each function
vspacing = 0.6; % vertical spacing between functions
vsize = 0.4;    % vertical size
vgap = vspacing - vsize;

for d = 0:dmax
  c = vsize*FindP(d, j);
  [m1, m2] = size(c);

  clear f;
  for k = 1:m2
    f(k,:) = EvalCombo(d, j, c(:,k), n) + vgap/2 + vspacing*(m2/2 - k);
  end;

  x = hsize*XVec(d, j, n) + d*hspacing;
  xl = [0, 0, hsize] + d*hspacing;
  yl = vgap/2 + vspacing*(m2/2 - [[1:m2]-vsize/vspacing; 1:m2; 1:m2]);
  plot(xl, yl, 'w');
  plot(x, f, 'w', 'LineWidth', 2);
end;

hold off;
axis equal;
axis off;

if nargin > 2
  filename = sprintf('spline-all-scaling-lev%i.ps', j-1);
  disp(['Saving EPS to "' filename '"...']);
  print('-deps', filename);
  % Latex stuff:
  disp(['    \begin{picture}(' num2str(dmax*hspacing + hsize) ',' ...
	  num2str(vspacing*(m2 + 0.5) - vgap) ')(0,' ...
	  num2str(vgap/2 - vspacing*(m2/2 + 0.5)) ')']);
  disp(['      \putbl{0,' num2str(vgap/2 - vspacing*m2/2) ...
	  '}{\psfig{file=' filename ',width=' ...
	  num2str(dmax*hspacing + hsize) 'in}}']);
  for d = 0:dmax
    m2 = 2^(j-1)+d;
    for k = 1:m2
      disp(['      \putr{' num2str(d*hspacing) ',' ...
	      num2str((vgap + vsize)/2 + vspacing*(m2/2 - k)) ...
	      '}{$N^' num2str(d) '_' num2str(k - 1) '\;$}']);
    end;
    disp(['      \putb{' num2str(d*hspacing + hsize/2) ',' ...
	    num2str(vgap/2 - vspacing*(m2/2 + 0.5)) ...
	    '}{\small degree ' num2str(d) '}']);
  end;
end;
return;

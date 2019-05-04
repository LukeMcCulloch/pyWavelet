function HaarApprox(normalization, save)

% HaarApprox(normalization, save) does a Haar wavelet decomposition of a
% simple function, and plots the scaling function and wavelet coefficients.
% The normalization parameter can be 'max' or 'L2'.  If save is present, the
% plot is saved as an EPS file.

% Settings:
filename = 'haar1d-approx.ps';
arrows = 0;
jinf = 8;
jmax = 4;
jmin = 0;
hspacing = 2.0; % horizontal spacing between functions
hsize = 1.6;    % horizontal size of each function
vspacing = 1.2; % vertical spacing between functions
vsize = 0.5;    % vertical size
vgap = vspacing - vsize;

% Set up plot:
clg;
hold on;
axis off;
axis equal;

% Axes:
axis1_x = [0, 0, 1];
axis1_y = [1, 0, 0];
axis2_x = [0, 0, nan, 0, 1];
axis2_y = [0.5, -0.5, nan, 0, 0];

% True function:
xinf = [0:2^jinf] / 2^jinf;
yinf = (10 * exp(-16*xinf/3) .* sin(8*xinf) + 1)/6;
plot(axis1_x*hsize, axis1_y*vsize + vspacing*(jmax-jmin), 'w');
plot(xinf*hsize, yinf*vsize + vspacing*(jmax-jmin), ':w');

% Finest scaling function approximation, computed by averaging:
y = zeros(1, 2^jmax);
for i = 1:2^(jinf-jmax)
  y = y + yinf(i:2^(jinf-jmax):2^jmax*(2^(jinf-jmax)-1)+i);
end;
y = y / 2^(jinf-jmax);
if strcmp(normalization, 'L2')
  y = y/2^(jmax/2);
end;
[xp, yp] = HaarCoords(y, normalization);
plot(xp*hsize, yp*vsize + vspacing*(jmax-jmin), 'w', 'LineWidth', 2);

if nargin > 1
  % Latex stuff:
  disp(['    \begin{picture}(' num2str(hspacing + hsize) ',' ...
	  num2str(vspacing*(jmax-jmin+1) - vgap/2) ')(0,' ...
	  num2str(-vgap/2) ')']);
  disp(['      \putbl{0,0}{\psfig{file=' filename ',width=' ...
	  num2str(hspacing + hsize) 'in}}']);
  disp(['      \putt{' num2str(hsize/2) ',' ...
	  num2str(vspacing*(jmax-jmin-0.1)) ...
	  '}{\small $V^' num2str(jmax) '$ approximation}']);
else
  plot(hsize/2, vspacing*(jmax-jmin)-0.1, 'o');
end;

% Repeatedly do summing and differencing on shorter and shorter part of y.
j = jmax;
len = 2^j;
while len >= 2
  y = HaarDecompStep(y, len, normalization);
  len = len/2;
  j = j - 1;
  [xp, yp] = HaarCoords(y(1:len), normalization);
  plot(axis1_x*hsize, axis1_y*vsize + vspacing*(j-jmin), 'w');
  plot(xinf*hsize, yinf*vsize + vspacing*(j-jmin), ':w');
  plot(xp*hsize, yp*vsize + vspacing*(j-jmin), 'w', 'LineWidth', 2);
  if nargin > 1
    % Latex stuff:
    if arrows
      disp(['      \put(' num2str(hsize/2) ',' ...
	      num2str(vspacing*(j-jmin+0.75)) ...
	      '){\vector(0,-1){' num2str(0.7*vspacing-vsize) '}}']);
    end;
    disp(['      \putt{' num2str(hsize/2) ',' ...
	    num2str(vspacing*(j-jmin)-0.1) ...
	    '}{\small $V^' num2str(j) '$ approximation}']);
  else
    if arrows
      plot([hsize/2 hsize/2], ...
	  [vspacing*(j-jmin+0.75) ...
	      vspacing*(j-jmin+0.75)-(0.7*vspacing-vsize)]);
    end;
    plot(hsize/2, vspacing*(j-jmin)-0.1, 'o');
  end;
  xp = [1:2:2^(j+1)-1]/2^(j+1);
  yp = y(len+1:2*len);
  plot(axis2_x*hsize + hspacing, ...
      axis2_y*vsize + vspacing*(j-jmin) + vsize/2, 'w');
  plot(xp*hsize + hspacing, 2*yp*vsize + vspacing*(j-jmin) + vsize/2, 'wo');
  w = 0.5*hsize;
  h = 0.7*vspacing-vsize;
  slope = round(w/h);
  w = h*slope;
  if nargin > 1
    % Latex stuff:
    if arrows
      disp(['      \put(' num2str(0.9*hsize) ',' ...
	      num2str(vspacing*(j-jmin+0.75)) ...
	      '){\vector(' num2str(slope) ...
	      ',-1){' num2str(w) '}}']);
    end;
    disp(['      \putt{' num2str(hsize/2 + hspacing) ',' ...
	    num2str(vspacing*(j-jmin)-0.1) ...
	    '}{\small $W^' num2str(j) ...
	    '$ detail coefficient' (j>0)*'s}'+(j==0)*'} ']);
  else
    if arrows
      plot([0.9*hsize 0.9*hsize+w], ...
	  [vspacing*(j-jmin+0.75) vspacing*(j-jmin+0.75)-(w/slope)]);
      end;
    plot(hsize/2 + hspacing, vspacing*(j-jmin)-0.1, 'o');
  end;
end;

if nargin > 1
  disp(['Saving EPS to "' filename '"...']);
  print('-deps', filename);
end;

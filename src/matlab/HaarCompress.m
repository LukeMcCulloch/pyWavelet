function HaarCompress(normalization, save)

% HaarCompress(normalization, save) does a Haar wavelet decomposition of a
% simple function, sorts the coefficients, and plots compression attempts.
% The normalization parameter can be 'max' or 'L2'.  If save is present, the
% plot is saved as an EPS file.

% Settings:
filename = 'haar1d-compress.ps';
jinf = 8;
jmax = 4;
jmin = 0;
hspacing = 2.0; % horizontal spacing between functions
hsize = 1.6;    % horizontal size of each function
vspacing = 0.9; % vertical spacing between functions
vsize = 0.5;    % vertical size
vgap = vspacing - vsize;

% Set up plot:
figure(1);
clg;
hold on;
axis off;
axis equal;

% Axis:
axis1_x = [0, 0, 1];
axis1_y = [1, 0, 0];

% True function:
xinf = [0:2^jinf] / 2^jinf;
yinf = (10 * exp(-16*xinf/3) .* sin(8*xinf) + 1)/6;

% Finest scaling function approximation, computed by averaging:
y = zeros(1, 2^jmax);
for i = 1:2^(jinf-jmax)
  y = y + yinf(i:2^(jinf-jmax):2^jmax*(2^(jinf-jmax)-1)+i);
end;
y = y / 2^(jinf-jmax);

% Set up plot:
cmax = 2^jmax;
cmin = 2;
cstep = 2;
rows = (cmax - cmin + cstep) / (2*cstep);
clg;
hold on;
axis off;
axis equal;

% Sort transformed coefficients:
xform = HaarDecomposition(y, normalization);
[junk, ind] = sort(abs(xform));

% Latex stuff:
if nargin > 1
  % Latex stuff:
  disp(['    \begin{picture}(' num2str(hspacing + hsize) ',' ...
	  num2str(vspacing*rows - vgap/2) ')(0,' num2str(-vgap/2) ')']);
  disp(['      \putbl{0,0}{\psfig{file=' filename ',width=' ...
	  num2str(hspacing + hsize) 'in}}']);
end;

% Show a bunch of approximations:
s = 0;
for c = cmax:-cstep:cmin
  % use only c coefficients
  xform(ind(1:2^jmax-c)) = zeros(1,2^jmax-c);
  y = HaarReconstruction(xform, normalization);
  [xp, yp] = HaarCoords(y, 'max');
  plot(axis1_x*hsize + rem(s,2)*hspacing, ...
      axis1_y*vsize + (rows - fix(s/2) - 1)*vspacing, 'w');
  plot(xinf*hsize + rem(s,2)*hspacing, ...
      yinf*vsize + (rows - fix(s/2) - 1)*vspacing, ':w');
  plot(xp*hsize + rem(s,2)*hspacing, ...
      yp*vsize + (rows - fix(s/2) - 1)*vspacing, 'w', 'LineWidth', 2);
  if nargin > 1
    % Latex stuff:
    disp(['      \putt{' num2str(hsize/2 + rem(s,2)*hspacing) ',' ...
	    num2str(vspacing*(rows - fix(s/2) - 1.1)) ...
	    '}{\small ' num2str(c) ' out of ' num2str(2^jmax) ...
	    ' coefficients}']);
  end;
  s = s + 1;
end;

if nargin > 1
  disp(['Saving EPS to "' filename '"...']);
  print('-deps', filename);
end;

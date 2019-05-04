function PlotAll(d, j, save)

% PlotAll(d, j, save) plots B-spline scaling functions and wavelets of
% degree d, level j.  If save is present, plot is saved as EPS file.

n = 20;         % samples per interval
hspacing = 2.5; % horizontal spacing between functions
hsize = 2;      % horizontal size of each function
vspacing = 0.4; % vertical spacing between functions

pscale = [0.4 0.7 0.7 0.7]*vspacing;
qscale = [0.4 0.5 0.9 1.1]*vspacing;

P = FindP(d, j+1);
P = P/max(max(abs(P)))*pscale(d+1);
Q = FindQ(d, j+1, 'L2');
Q = Q/max(max(abs(Q)))*qscale(d+1);

figure(1);
clg;
hold on;

[m1, m2] = size(P);
x = hsize*XVec(d, j+1, n);
for k = 1:m2
  f(k,:) = EvalCombo(d, j+1, P(:,k), n) + vspacing*(m2 - k);
end;
plot([0, hsize], vspacing*[0:m2-1; 0:m2-1], 'w');
plot(x, f, 'w', 'LineWidth', 2);
num = m2;

[m1, m2] = size(Q);
x = x + hspacing;
clear f;
for k = 1:m2
  f(k,:) = EvalCombo(d, j+1, Q(:,k), n)+ vspacing*((num + m2)/2 - k);
end;
plot(hspacing + [0, hsize], vspacing*([0:m2-1; 0:m2-1] + (num - m2)/2), 'w');
plot(x, f, 'w', 'LineWidth', 2);

hold off;
axis equal;
axis off;

if nargin > 2
  filename = sprintf('spline-deg%i-lev%i.ps', d, j);
  disp(['Saving EPS to "' filename '"...']);
  print('-deps', filename);
end;
return;

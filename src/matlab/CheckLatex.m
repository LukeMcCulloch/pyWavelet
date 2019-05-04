% CheckLatex checks the Qs printed by LatexQ for normalization and
% orthogonality to P'*I.

for d = 0:3
  for j = 1:4
    disp(['---------- d = ' num2str(d) ', j = ' num2str(j) ' -----------']);
    P = FindP(d, j);
    P_denom = LCD(P);
    P = round(P*P_denom);
    disp(['P = 1/' num2str(P_denom)]);
    disp(P);
    I = Inner(d, j);
    I_denom = LCD(I);
    I = round(I*I_denom);
    disp(['I = 1/' num2str(I_denom)]);
    disp(I);
    disp('Q = ');
    Q = LatexQ(d, j)
    ip = P'*I*Q/I_denom;
    if all(all(abs(ip) < 1e-5))
      disp('******* Passed orthogonality test.');
    else
      disp(['******* FAILED orthogonality: IP < ' num2str(max(max(abs(ip))))]);
    end;
    ip = diag(Q'*I*Q/I_denom);
    if all(abs(ip - 1) < 1e-5)
      disp('******* Passed normalization test.');
    else
      disp(['******* FAILED normalization: |IP-1| < ' ...
	      num2str(max(abs(ip - 1)))]);
    end;
    disp('Press a key...');
    pause;
  end;
end;
return;

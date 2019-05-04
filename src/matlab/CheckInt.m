% CheckInt checks integer versions of Q's for orthogonality to P'*I.

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
    Q = FindQ(d, j, 'lcd')
    ip = P'*I*Q;
    if all(all(ip == 0))
      disp('******* Passed orthogonality test.');
    else
      disp('******* FAILED!');
    end;
    disp('Press a key...');
    pause;
  end;
end;
return;

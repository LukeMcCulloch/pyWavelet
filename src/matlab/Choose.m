function c = Choose(n, r)

% c = Choose(n, r) returns (n choose r) = n! / (r! (n-r)!).

c = Factorial(n)./(Factorial(r).*Factorial(n-r));
return;

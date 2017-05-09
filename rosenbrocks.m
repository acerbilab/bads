function y = rosenbrocks(x)
%ROSENBROCKS Rosenbrock's 'banana' function in any dimension.

y = sum(100*(x(:,1:end-1) .^2 - x(:,2:end)) .^ 2 + (x(:,1:end-1)-1) .^2);
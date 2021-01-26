function [x, success] = nonlinear_picard2(getf, x, tol, kmax)
% getf = anonymous function, defined as e.g. getf = @(x) x^2 - 3*x + 2
% x = scalar, initial guess for the root
% tol = scalar, criterion for convergence
% kmax = maximum number of iterations

f = getf(x); %evaluate getf at initial guess
k = 0; %initialize counter to 0

while abs(f) > tol
   x = x + f(x); %Picard's update rule
   k = k + 1; %update counter by +1
   f = getf(x); % update function value
   if k > kmax || x == Inf || x == -Inf
      warning('Did not converge \n')
      break
   end
end

if k > kmax || x == Inf || x == -Inf
    success = 0;
else
    success = 1;
end
end

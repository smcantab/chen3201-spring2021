clear;clc

getf = @(x) sin(pi*x);
getdf = @(x) pi * cos(pi*x);

xinit = linspace(0, 2, 201);
roots = zeros(length(xinit), 1);
success = zeros(length(xinit), 1);
tol = 1e-10;
maxcount = 10000;

for i = 1:length(xinit)
   [roots(i), success(i)] = nonlinear_newton(getf, getdf, xinit(i), tol, maxcount);
   if ~success(i) %check if successful, we expect the method to fail when rem(x, 0.5) = 0
        roots(i) = [];
        success(i) = [];
        fprintf('deleted value %.3f', roots(i));
   end
end

plot(xinit, roots, '-o', 'MarkerSize', 10)
hold on
plot(xinit, (arrayfun(getf, xinit)), '-r', xinit, (arrayfun(getdf, xinit)), '--b', 'LineWidth', 2)
hold off
grid on, grid minor
xlabel('Initial guess for x'), ylabel('Root from Newton Method')
legend('roots', '|f|', '|df|')
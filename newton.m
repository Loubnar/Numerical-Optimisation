function [x, iter, temps, exiter] = newton(f,grad_f,hess_f,x0,epsilon,i)
%% title : Local Newton Method
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% f       : Function to evaluate
% grad_f  : Function's gradient
% hess_f  : Function's hessian
% x0      : Initial value of x
% epsilon : Criteria precision
%% Output
% x       : Solution's approximation
% iter    : Number of iterations
% time    : Computation time
% i       : Integer indicating the method of search
%% 
global nfev ngev nhev;
nfev=0;
ngev=0;
nhev=0;
tic;
x= x0;
x_old=x0;
iter=0;
s = 1;
rho = 0.9;
exiter = 0;
temp = 0;
while (norm(grad_f(x)+hess_f(x)*(x-x_old))> epsilon) || (iter==0)
    % While the gradient of the quadratic is not close to zero...
    d = -hess_f(x)\grad_f(x);
	switch i
	case 1
		s=1;
	case 2
		[s, vals, valfs, temp] = backtracking(f, grad_f, x, d, rho);
	case 3
		[s, vals, valfs, temp] = bissection(f, grad_f, x, d);
	case 4
		[s, vals, valfs, temp] = interpolation(f, grad_f, x, d, s);
	case 5
		[s, vals, valfs, temp] = approach(f, grad_f, x, d, s, 10);
	end
    x_old = x;
    x = x + s*d;
    iter = iter +1;
	exiter = exiter+temp;
end
time = toc;

    
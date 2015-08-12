function [x, iter, time, exiter] = bfgs(f,grad_f,hess_f,x0,epsilon,i)
%% title : Quasi-Newton BFGS method
%% author : Loubna Rizqi
%% date : 21/01/2015
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
%% Implantation de la methode de BFGS
tic;
global nfev ngev nhev;
nfev=0;
ngev=0;
nhev=0;
x = x0;
iter=0;
x_old = x0;
H = inv(hess_f(x));
s = 1;
exiter = 0;
temp = 0;
rho=0.9;
while ((norm(grad_f(x))> epsilon) || (iter==0))
    d = -H * grad_f(x);
    x_old = x;
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
    x = x + s*d;
    y = grad_f(x) - grad_f(x_old);
    H = (eye(size(d,1))-(d*y')/(d'*y))*H*(eye(size(d,1))-(y*d')/(d'*y))+s*(d*d')/(d'*y);
    iter = iter +1;
	exiter = exiter +temp;
end
time=toc;
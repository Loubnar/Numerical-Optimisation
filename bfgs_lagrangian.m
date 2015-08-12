function [x, iter] = bfgs_lagrangian(f,grad_f,x0,epsilon)
%% title : Quasi-Newton BFGS method with backtracking for the augmented lagrangian algorithm
%% author : Loubna Rizqi
%% date : 21/01/2015
%% Input
% f       : Function to evaluate
% grad_f  : Function's gradient
% hess_f  : Function's hessian
% x0      : Initial value of x
% epsilon : Criteria precision
% c       : Constraints
% J_c     : Constraints jacobian
%% Output
% x       : Solution's approximation
% iter    : Number of iterations
% time    : Computation time
% exiter  : Number of external iterations
%% Implementation of BFGS with Backtracking
global nfev ngev nhev;
x = x0;
iter=0;
x_old = x0;
H = eye(2);
s = 1;
rho = 0.9;
exiter=0;
while (norm(grad_f(x))> epsilon)
    d = -H * grad_f(x);
    [s, vals, valfs, temp] = backtracking(f, grad_f, x, d, rho);
    x_old = x;
    x = x + s*d;
    y = grad_f(x) - grad_f(x_old);
    H = (eye(size(d,1))-(d*y')/(d'*y))*H*(eye(size(d,1))-(y*d')/(d'*y))+s*(d*d')/(d'*y);
    iter = iter +1;
    exiter = exiter + temp;
end
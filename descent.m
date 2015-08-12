function [x, iter, time] = descent(grad_f,hess_f,x0,epsilon)
%% title : Descent Algorithm
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
%%
global nfev ngev nhev;
nfev=0;
ngev=0;
nhev=0;
tic;
x=x0;
iter=0;
s = 1; % Descent step
while (norm(grad_f(x)) > epsilon) || (iter==0) % While the exit criteria is not met...
    d = -grad_f(x); %... recalculate the steepest descent direction
    s = -(d'*grad_f(x))/(d'*hess_f(x)*d); %... find a step which satisfies f(xk+1)<=f(xk)
    x = x + s * d; %... update x
    iter  = iter + 1;
end
time = toc;

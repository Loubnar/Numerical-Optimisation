function [x, iter, time, exiter] = lagrangian(f,grad_f, hess_f, c, J_c, x0, lambda0, tho, mu0)
%% title : Augmented Lagrangian Algorithm
%% author : Loubna Rizqi
%% date : 21/01/2015
%% Input 
% f       : Function to evaluate
% grad_f  : Function's gradient
% hess_f  : Function's hessian
% x0      : Initial value of x
% epsilon : Criteria precision
% c       : Constraints
% J_c     : Constraints' jacobian
%% Output 
% x       : Solution's approximation
% iter    : Number of iterations
% time    : Computation time
% exiter  : Number of external iterations
%%
tic;
global nfev ngev nhev nfcev ngcev nhcev;
nfev=0;
ngev=0;
nhev=0;
nfcev=0;
ngcev=0;
nhcev=0;
x = x0;
lambda = lambda0;
eta0 = 0.1258925;
alpha = 0.1;
beta = 0.9;
mu = mu0;
epsilon0 = 1/mu0;
epsilon = epsilon0;
etak = eta0/(mu0^alpha);
iter=0;
x = x0;
exiter = 0;
while true
         %%%% Searching a minimizer using bfgs %%%%
         [x, temp] = bfgs_lagrangian(@(x) f(x)+lambda'*c(x)+(mu/2)*(norm(c(x)))^2,@(x) grad_f(x)+lambda'*J_c(x)+mu*J_c(x)*c(x),x,epsilon);
         %%%% test sur la convergence %%%%
         if (norm(grad_f(x)+lambda'*J_c(x)+mu*J_c(x)*c(x)) <= 10^-6) && (norm(c(x))<=10^-6)
             break;
         end
         %%%% Update test %%%%
         if norm(c(x))<=etak % If the norm of the constraints is smaller than etak...
            %... Update the lagrangian multipliers
            lambda = lambda + mu*c(x);
            epsilon = epsilon/mu;
            etak = etak/(mu^beta);
         else 
            %...otherwise update paremeters of penalty
            mu = tho*mu;
            epsilon = epsilon0/mu;
            etak = eta0/(mu^alpha);
         end
         iter = iter+1;
         exiter = exiter +temp;
end
time =toc;
            
            
         

function [s, vals, valsf, iter] = backtracking(f, grad_f, x, d, rho, c1)
%% title : Backtracking algorithm
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% f       : Function to evaluate
% grad_f  : Function's gradient
% x       : Starting point
% d       : Descent direction
% rho     : Contraction direction
% c1      : coefficient of the first condition of Wolfe
%% Output 
% s       : Solution's approximation
% iter    : Number of iterations
% vals    : Step value
% valsf   : Values of f applied to x+s*d
%% 
global nfev ngev nhev;
s = 1;
vals = [s];
valsf = [f(x+d*s)];
iter = 0;
c1 = 10^-3;
while (f(x+s*d) > (f(x) + s*c1*grad_f(x)'*d)) && iter<50
    s = rho * s;
    iter = iter + 1;
    vals = [vals s];
    valsf = [valsf f(x+s*d)];
end    
